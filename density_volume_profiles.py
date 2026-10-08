# protein_volume_multipleconf_faster.py

import numpy as np
import pandas as pd
from collections import defaultdict
import os
from pathlib import Path
import time
import glob

# Define van der Waals radii for common elements 
van_der_waals_radii = {
    'H': 1.2,  # Example values; replace with accurate ones
    'O': 1.52,
    'C': 1.7,
    'N': 1.55,
    'P': 1.8,
}

list_atom_names = ['C', 'H', 'N', 'O', 'P']

def read_file(filename):
    atoms = []
    with open(filename, 'r') as file:
        for line in file:
          parts = line.split()
          if len(parts) > 3: 
            atom_name = parts[0]
            x, y, z = map(float, parts[1:])
            radius = van_der_waals_radii.get(atom_name, 1.5)  # Default radius if not specified
            atoms.append((atom_name, x, y, z, radius))
    return atoms

def read_xyz_trajectory(filename):
    """
    Read a multi-frame XYZ trajectory.

    Returns
    -------
    frames : list
        [
            [(atom_name, x, y, z, radius), ...],
            [(atom_name, x, y, z, radius), ...],
            ...
        ]
    """

    frames = []

    with open(filename, "r") as f:
        while True:
            line = f.readline()

            if not line:
                break

            line = line.strip()
            if not line:
                continue

            try:
                natoms = int(line)
            except ValueError:
                raise ValueError(f"Expected atom count line in XYZ trajectory, got: {line}")

            # comment line
            f.readline()

            atoms = []

            for _ in range(natoms):
                atom_line = f.readline()
                if not atom_line:
                    raise ValueError("Unexpected end of file while reading XYZ trajectory")

                parts = atom_line.split()
                if len(parts) < 4:
                    continue

                atom_name = parts[0]
                x = float(parts[1])
                y = float(parts[2])
                z = float(parts[3])
                radius = van_der_waals_radii.get(atom_name, 1.5)

                atoms.append((atom_name, x, y, z, radius))

            frames.append(atoms)

    return frames

def calculate_min_max(atoms):
    x_coords = [atom[1] for atom in atoms]
    y_coords = [atom[2] for atom in atoms]
    z_coords = [atom[3] for atom in atoms]

    min_x, max_x = min(x_coords), max(x_coords)
    min_y, max_y = min(y_coords), max(y_coords)
    min_z, max_z = min(z_coords), max(z_coords)
    return (min_x, max_x), (min_y, max_y), (min_z, max_z)

def count_cubes(atoms, min_max_coords, bin_size=0.5):
    (min_x, max_x), (min_y, max_y), (min_z, max_z) = min_max_coords

    # Calculate number of bins on each axis
    bins_x = int(np.ceil((max_x - min_x) / bin_size))
    bins_y = int(np.ceil((max_y - min_y) / bin_size))
    bins_z = int(np.ceil((max_z - min_z) / bin_size))

    # Initialize 3D grid for bins and a list for slices
    grid = np.zeros((bins_x, bins_y, bins_z), dtype=bool)
    slices = [defaultdict(set) for _ in range(bins_z)]
    slices_counts = [defaultdict(int) for _ in range(bins_z)]
    slice_occupied_cubes = [0] * bins_z
 
    for atom_name, x, y, z, radius in atoms:
        # Calculate the range of bins affected by the atom's van der Waals radius
        min_i = int(np.floor((x - radius - min_x) / bin_size))
        max_i = int(np.ceil((x + radius - min_x) / bin_size))
        min_j = int(np.floor((y - radius - min_y) / bin_size))
        max_j = int(np.ceil((y + radius - min_y) / bin_size))

   #    not applying the radius on z to avoid ending up into other slices
   #     min_k = int(np.floor((z - radius - min_z) / bin_size))
   #     max_k = int(np.ceil((z + radius - min_z) / bin_size))     
        min_k = int(np.floor((z - bin_size/2 - min_z) / bin_size))
        max_k = int(np.ceil((z + bin_size/2 - min_z) / bin_size))


        # Mark all affected bins and count cubes in each slice: *** with volumes ***
        for i in range(min_i, max_i):
            for j in range(min_j, max_j):
                for k in range(min_k, max_k):
                    if 0 <= i < bins_x and 0 <= j < bins_y and 0 <= k < bins_z:
                        if not grid[i, j, k]:  # Only count if the cube is not already occupied
                            grid[i, j, k] = True
                            slice_occupied_cubes[k] += 1
                        slices[k][atom_name].add((i, j, k))

        # Count each atom exactly once, including points on grid planes.
        # The last bin includes the upper box boundary.
        coords = (x, y, z)
        bounds = min_max_coords
        if all(lo <= value <= hi for value, (lo, hi) in zip(coords, bounds)):
            k = min(int(np.floor((z - min_z) / bin_size)), bins_z - 1)
            slices_counts[k][atom_name] += 1

    for counts in slices_counts:
        for element in list_atom_names:
            counts.setdefault(element, 0)
    return slice_occupied_cubes, slices_counts, bins_z


def calculate_volume_ND_from_atoms(atoms, box_size, bin_size, bounds=None):
    """
    Calculate volume fraction and number density profiles from an atom list.

    atoms format:
        [(atom_name, x, y, z, radius), ...]
    """

    actual_min_max_coords = calculate_min_max(atoms)

    if bounds is not None:
        min_max_coords = bounds
    elif box_size == None:
        min_max_coords = actual_min_max_coords
    else:
        min_max_coords = (
            (-box_size / 2, box_size / 2),
            (-box_size / 2, box_size / 2),
            (-box_size / 2, box_size / 2)
        )

    if actual_min_max_coords[0][0] < min_max_coords[0][0] or actual_min_max_coords[0][1] > min_max_coords[0][1]:
        print("WARNING: Box %s %s not large enough along X axis!" % (min_max_coords[0][0], min_max_coords[0][1]))
    if actual_min_max_coords[1][0] < min_max_coords[1][0] or actual_min_max_coords[1][1] > min_max_coords[1][1]:
        print("WARNING: Box %s %s not large enough along Y axis!" % (min_max_coords[1][0], min_max_coords[1][1]))
    if actual_min_max_coords[2][0] < min_max_coords[2][0] or actual_min_max_coords[2][1] > min_max_coords[2][1]:
        print("WARNING: Box %s %s not large enough along Z axis!" % (min_max_coords[2][0], min_max_coords[2][1]))

    area = (min_max_coords[0][1] - min_max_coords[0][0]) * (min_max_coords[1][1] - min_max_coords[1][0])
    slice_volume = area * bin_size

    # count cubes for volume calc
    slice_occupied_cubes, slices_counts, num_slices = count_cubes(atoms, min_max_coords, bin_size)

    total_volume = 0

    dict_vol_frac = {}
    dict_ND_C = {}
    dict_ND_H = {}
    dict_ND_N = {}
    dict_ND_O = {}

    for k in range(num_slices):
        occupied_volume_slice = slice_occupied_cubes[k] * bin_size * bin_size * bin_size
        total_volume = total_volume + occupied_volume_slice

        z = min_max_coords[2][0] + k * bin_size
        dict_vol_frac[z] = occupied_volume_slice / slice_volume

        for atom_type, count in slices_counts[k].items():
            if atom_type == 'C':
                dict_ND_C[z] = count / slice_volume
            elif atom_type == 'H':
                dict_ND_H[z] = count / slice_volume
            elif atom_type == 'N':
                dict_ND_N[z] = count / slice_volume
            elif atom_type == 'O':
                dict_ND_O[z] = count / slice_volume

    return dict_vol_frac, total_volume, dict_ND_C, dict_ND_H, dict_ND_N, dict_ND_O


def collect_profile_values(values):
    """
    Convert one calculate_volume_ND result into simple profile lists.
    """

    volume_fraction = values[0]
    volume = values[1]
    density_C = values[2]
    density_H = values[3]
    density_N = values[4]
    density_O = values[5]

    keys = []
    C_values = []
    for key, val in density_C.items():
        keys.append(key)
        C_values.append(val)

    H_values = []
    for key, val in density_H.items():
        H_values.append(val)

    N_values = []
    for key, val in density_N.items():
        N_values.append(val)

    O_values = []
    for key, val in density_O.items():
        O_values.append(val)

    VF_values = []
    for key, val in volume_fraction.items():
        VF_values.append(val)

    return keys, C_values, H_values, N_values, O_values, VF_values, volume

def average_profiles(carbon, hydrogen, nitrogen, oxygen, vol_frac, volumes_list):
    """
    Average profiles over conformations.

    Rows = conformations
    Columns = z-bins
    """

    c = pd.DataFrame(carbon)
    h = pd.DataFrame(hydrogen)
    n = pd.DataFrame(nitrogen)
    o = pd.DataFrame(oxygen)
    vf = pd.DataFrame(vol_frac)

    averaged = {
        "C_mean": c.mean(axis=0),
        "C_std": c.std(axis=0).fillna(0.0),
        "H_mean": h.mean(axis=0),
        "H_std": h.std(axis=0).fillna(0.0),
        "N_mean": n.mean(axis=0),
        "N_std": n.std(axis=0).fillna(0.0),
        "O_mean": o.mean(axis=0),
        "O_std": o.std(axis=0).fillna(0.0),
        "VF_mean": vf.mean(axis=0),
        "VF_std": vf.std(axis=0).fillna(0.0),
        "average_volume": sum(volumes_list) / len(volumes_list),
    }

    return averaged

def write_profile_file(filename, keys, mean_values, std_values, open_flag="w"):
    """Write z, mean, std profile file."""
    with open(filename, open_flag) as g:
        for i in range(0, len(keys)):
            g.write(str(keys[i]) + " " + str(mean_values[i]) + " " + str(std_values[i]) + "\n")

def calculate_volume_ND(filename, box_size, bin_size):
    atoms = read_file(filename)
    return calculate_volume_ND_from_atoms(atoms, box_size, bin_size)


def calcAll(size,w_slice, hydrogen_group2=None):
    start_time = time.time()
    split_profiles = {"H_group1": [], "H_group2": []}
    carbon = []
    hydrogen = []
    oxygen = []
    nitrogen = []
    vol_frac = []
    volumes_list = []
    sourcedir = Path("tempFiles")

    # Number densities written to density_*.dat. Do we append?
    append_output = False
    if append_output:
        open_flag = 'a'
    else:
        open_flag = 'w'

    trajectory_file = sourcedir / "trajectory.xyz"

    if trajectory_file.exists():
        # New mode:
        # one multi-frame XYZ, one frame per conformation
        frames = read_xyz_trajectory(trajectory_file)

        inputs = []
        use_frames = True
        for frame in frames:
            inputs.append(frame)


    else:
        # Old mode:
        # one XYZ file per conformation
        inputs = glob.glob(str(sourcedir) + '/*.xyz')
        use_frames = False


    if not inputs:
        raise ValueError("No coordinate conformations found in tempFiles.")
    if not use_frames:
        inputs = [read_file(filename) for filename in sorted(inputs)]
    if not np.isfinite(w_slice) or w_slice <= 0:
        raise ValueError("Bin size must be finite and positive.")
    if size is not None and (not np.isfinite(size) or size <= 0):
        raise ValueError("Box size must be finite and positive.")
    if hydrogen_group2 is not None and len(hydrogen_group2) != len(inputs):
        raise ValueError("Expected one hydrogen mask per conformation.")
    # A common grid is essential when averaging conformations without a box.
    if size is None:
        coords = np.concatenate([np.array([a[1:4] for a in frame]) for frame in inputs])
        lo = np.floor(coords.min(axis=0) / w_slice) * w_slice
        hi = np.maximum(np.ceil(coords.max(axis=0) / w_slice) * w_slice, lo + w_slice)
        bounds = tuple(zip(lo, hi))
    else:
        bounds = ((-size / 2, size / 2),) * 3

    # loop over conformations
    # In new mode, each item is a frame from trajectory.xyz.
    # In old mode, each item is one XYZ filename.
    n_inputs = len(inputs)

    for input_index, item in enumerate(inputs, start=1):
        print(f"[pvm] processing conformation {input_index}/{n_inputs}")

        values = calculate_volume_ND_from_atoms(item, box_size=size, bin_size=w_slice,
                                                bounds=bounds)
        if hydrogen_group2 is not None:
            mask = np.asarray(hydrogen_group2[input_index - 1], dtype=bool)
            is_h = np.array([a[0] == 'H' for a in item])
            if mask.shape != is_h.shape or np.any(mask & ~is_h):
                raise ValueError("Hydrogen mask must match atom order and select only H.")
            n_bins = len(values[3])
            slice_volume = (bounds[0][1]-bounds[0][0]) * (bounds[1][1]-bounds[1][0]) * w_slice
            for name, selected in (("H_group1", is_h & ~mask), ("H_group2", mask)):
                counts = np.zeros(n_bins)
                for index in np.flatnonzero(selected):
                    xyz = item[index][1:4]
                    if all(lo <= v <= hi for v, (lo, hi) in zip(xyz, bounds)):
                        k = min(int(np.floor((xyz[2]-bounds[2][0])/w_slice)), n_bins-1)
                        counts[k] += 1
                split_profiles[name].append(counts / slice_volume)

        volume_fraction = values[0]
        volume = values[1]
        density_C = values[2]
        density_H = values[3]
        density_N = values[4]
        density_O = values[5]

        C_values = []
        keys = []
        for key, val in density_C.items():
            C_values.append(val)
            keys.append(key)
        carbon.append(C_values)

        H_values = []
        for key, val in density_H.items():
            H_values.append(val)
        hydrogen.append(H_values)

        N_values = []
        for key, val in density_N.items():
            N_values.append(val)
        nitrogen.append(N_values)

        O_values = []
        for key, val in density_O.items():
            O_values.append(val)
        oxygen.append(O_values)

        volume_fraction_values = []
        for key, val in volume_fraction.items():
            volume_fraction_values.append(val)
        vol_frac.append(volume_fraction_values)

        volumes_list.append(volume)

    # Having collated ND profiles for each conformation, calculate mean profile.
    # This assumes consistent bins across the different input conformations.
    # Pandas dataframe has files as rows (axis=0) and bins as columns (axis=1).
    # Calculate mean and standard deviation across files, for each z-bin.

    averaged = average_profiles(
        carbon=carbon,
        hydrogen=hydrogen,
        nitrogen=nitrogen,
        oxygen=oxygen,
        vol_frac=vol_frac,
        volumes_list=volumes_list,
    )

    print("average volume=", averaged["average_volume"])

    write_profile_file("density_C.dat", keys, averaged["C_mean"], averaged["C_std"], open_flag)
    write_profile_file("density_H.dat", keys, averaged["H_mean"], averaged["H_std"], open_flag)
    write_profile_file("density_N.dat", keys, averaged["N_mean"], averaged["N_std"], open_flag)
    write_profile_file("density_O.dat", keys, averaged["O_mean"], averaged["O_std"], open_flag)
    write_profile_file("vol_frac_along_z.dat", keys, averaged["VF_mean"], averaged["VF_std"], open_flag)

    if hydrogen_group2 is not None:
        for name, profiles in split_profiles.items():
            frame = pd.DataFrame(profiles)
            write_profile_file(f"density_{name}.dat", keys, frame.mean(axis=0),
                               frame.std(axis=0).fillna(0.0), open_flag)
        if not np.allclose(np.array(split_profiles["H_group1"]) +
                           np.array(split_profiles["H_group2"]), np.array(hydrogen)):
            raise RuntimeError("Hydrogen partition does not conserve total density.")

    elapsed = time.time() - start_time
    print(f"[pvm] finished profile calculation in {elapsed:.2f} s")

def read_three_col(path):
    """Read z, mean, std (or ignore lines that don't parse)."""
    z = []; mu = []; sd = []
    try:
        with open(path) as f:
            for line in f:
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        z.append(float(parts[0]))
                        mu.append(float(parts[1]))
                        sd.append(float(parts[2]) if len(parts) >= 3 else 0.0)
                    except ValueError:
                        continue
    except FileNotFoundError:
        pass
    return z, mu, sd

def return_densities_and_volumes(input_dir, box_size=None, bin_size=0.5, hydrogen_group2=None):
    """
    Importable function (for multiconf_tilt_table). Runs the existing calcAll on 'input_dir'.

    Args:
        input_dir: orientation directory containing tempFiles/*.xyz
        box_size: passed to calcAll(size=...)
        bin_size: passed to calcAll(w_slice=...)

    Returns:
        A dictionary looking like
            {
              "z": z_grid,
              "density": {"C":..., "H":..., "N":..., "O":..., "P":... (if present)},
              "vf": volume_fraction_profile,
            }
    """
    input_dir = Path(input_dir).resolve()
    if not input_dir.exists():
        raise FileNotFoundError(input_dir, "not found")

    # Run calcAll in that directory so outputs are written there
    import os
    cwd0 = os.getcwd()
    try:
        os.chdir(str(input_dir))
        # call the module's existing function (do not rename/change your calcAll)
        calcAll(size=box_size, w_slice=bin_size, hydrogen_group2=hydrogen_group2)
    finally:
        os.chdir(cwd0) # go back


    # Read outputs and store in dictionary
    dens = {}
    zC, C, _ = read_three_col(input_dir / "density_C.dat") # get also the z once for all
    dens["C"] = C
    _ , H, _ = read_three_col(input_dir / "density_H.dat")
    dens["H"] = H
    _ , N, _ = read_three_col(input_dir / "density_N.dat")
    dens["N"] = N
    _ , O, _ = read_three_col(input_dir / "density_O.dat")
    dens["O"] = O
    zVF, VF, _ = read_three_col(input_dir / "vol_frac_along_z.dat")

    # Optional phosphorus if written/present in system
    if (input_dir / "density_P.dat").exists():
        _ , P, _ = (read_three_col(input_dir / "density_P.dat"))
        dens["P"] = P

    if hydrogen_group2 is not None:
        for group in (1, 2):
            _, values, _ = read_three_col(input_dir / f"density_H_group{group}.dat")
            dens[f"H_group{group}"] = values

    return {"z": zC, "density": dens, "vf": VF}

# Argument parsing
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("input_dir", nargs="?", default=".")
    parser.add_argument("--box-size", type=float, default=None)
    parser.add_argument("--bin", type=float, default=0.5, dest="bin_size")
    args = parser.parse_args()
    return_densities_and_volumes(args.input_dir, box_size=args.box_size, bin_size=args.bin_size)

 
