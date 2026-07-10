#!/usr/bin/env python3
"""
multiconf_tilt_table_wrigglefirst.py — full pipeline for 1) Rotation 2) Wriggling 3) Calculation of ND 4) Calculation of SLD 5) Building a final array (lookup table)

Example:

python multiconf_tilt_table_wrigglefirst.py \
  --pdb peptide.pdb \
  --theta-angles 0:90:10 \
  --phi-angles 0:90:10 \
  --outdir orientations \
  --json my_lookup.json

Recent changes:
1) Option for lookup table custom name 
2) Automatic alignment and centering on origin of the input file before doing anything else
3) First rotation around theta, second rotation around phi (use of polar coordinates)

"""

import argparse, os, sys, subprocess
from pathlib import Path
import json
# Import from modules 
from io_utils import write_xyz_trajectory, pdb_to_atoms, write_pdb_from_template
import density_volume_profiles as dvp
from geometry_utils import center_and_align_atoms_along_x, center_atoms_on_origin, rotate_atoms_theta_phi
from density_utils import read_density_dat, compute_sld_profiles
from workflow_utils import clean_distributions, parse_angle_spec, ensure_stochastic_in_dir

# ------------------------ main API / pipeline helpers ------------------------

def write_run_metadata(out_root, pdb, theta_angles, phi_angles, box_size,
                       bin_size, nstruc, skip_wriggle, json_name, wriggle_mode):
    """Write a small metadata file describing the current run."""
    metadata = {
        "pdb": str(pdb),
        "theta_angles": str(theta_angles),
        "phi_angles": str(phi_angles),
        "box_size": box_size,
        "bin_size": bin_size,
        "nstruc": nstruc,
        "skip_wriggle": skip_wriggle,
        "json_name": json_name,
        "wriggle_mode": wriggle_mode,
    }

    metadata_path = out_root / "run_metadata.json"
    with open(metadata_path, "w") as f:
        json.dump(metadata, f, indent=2)

    return metadata_path

def prepare_wriggled_structures(
    pdb,
    out_root,
    python=sys.executable,
    skip_wriggle=False,
    wriggle_arg_input="--input",
    wriggle_arg_nstruc="-n",
    nstruc=20,
    wriggle_mode="protein",
    bb_sd=None,
    bb_window=None,
    bb_attempts=None,
    bb_stride=None,
    bb_clash_scale=None,
    sc_sd_backbone=None,
    sc_opt_passes=None,
    sc_opt_scale=None,
    sc_opt_var=None,

    ):
    """
    Prepare the aligned starting structure and generate wriggled conformations.

    This step is done only once, in a shared wriggles_base directory.
    All orientation folders later reuse the same wriggled structures.
    """

    # Wriggle only once, in a shared base directory
    wriggles_base = out_root / "wriggles_base"
    wriggles_base.mkdir(parents=True, exist_ok=True)

    src_pdb = Path(pdb)
    dst_pdb = wriggles_base / src_pdb.name

    # Center starting conformation on origin and align along x BEFORE any other step
    # Write this aligned structure into wriggles_base using the original filename,
    # so wriggle_driver.py (and downstream wriggle code) uses the aligned starting structure.
    atoms0 = pdb_to_atoms(src_pdb)
    atoms0_aligned = center_and_align_atoms_along_x(atoms0)
    write_pdb_from_template(src_pdb, dst_pdb, atoms0_aligned)

    # Make stochastic_wriggle.py available inside the wriggle directory
    ensure_stochastic_in_dir(wriggles_base)

    if skip_wriggle == False:
        wrigglefile_path = str((Path(__file__).parent / "wriggle_driver.py").resolve())

        wriggle_overrides = []

        if bb_sd is not None:
            wriggle_overrides += ["--bb-sd", str(bb_sd)]

        if bb_window is not None:
            wriggle_overrides += ["--bb-window", str(bb_window)]

        if bb_attempts is not None:
            wriggle_overrides += ["--bb-attempts", str(bb_attempts)]

        if bb_stride is not None:
            wriggle_overrides += ["--bb-stride", str(bb_stride)]

        if bb_clash_scale is not None:
            wriggle_overrides += ["--bb-clash-scale", str(bb_clash_scale)]

        if sc_sd_backbone is not None:
            wriggle_overrides += ["--sc-sd-backbone", str(sc_sd_backbone)]

        if sc_opt_passes is not None:
            wriggle_overrides += ["--sc-opt-passes", str(sc_opt_passes)]

        if sc_opt_scale is not None:
            wriggle_overrides += ["--sc-opt-scale", str(sc_opt_scale)]

        if sc_opt_var is not None:
            wriggle_overrides += ["--sc-opt-var", str(sc_opt_var)]


        subprocess.run(
            [
                python,
                str(wrigglefile_path),
                str(wriggle_arg_input),
                str(src_pdb.name),
                str(wriggle_arg_nstruc),
                str(nstruc),
                "--mode", str(wriggle_mode),
            ] + wriggle_overrides,
            cwd=str(wriggles_base),
            check=True,
        )

    wriggled_pdbs = sorted((wriggles_base / "structures").glob("*.pdb"))

    if len(wriggled_pdbs) == 0:
        raise RuntimeError(
            "No wriggled PDB files were found. "
            "Expected files under: "
            f"{wriggles_base / 'structures'}"
        )

    return wriggled_pdbs

def rotate_wriggled_structures_for_orientation(wriggled_pdbs, theta, phi, orient_dir):
    """
    Rotate all wriggled structures for one theta/phi orientation.

    The rotated conformations are written as one multi-frame XYZ trajectory:
        tempFiles/trajectory.xyz

    Each frame corresponds to one wriggled conformation.
    """

    temp_dir = orient_dir / "tempFiles"
    temp_dir.mkdir(exist_ok=True, parents=True)

    frames = []

    for pdbf in wriggled_pdbs:
        atoms = pdb_to_atoms(pdbf)

        # ensure each conformation is centered before rotations
        atoms = center_atoms_on_origin(atoms)

        # rotate: theta is tilt relative to z, phi is spin around the protein axis
        atoms_rot = rotate_atoms_theta_phi(atoms, theta, phi)

        # store this conformation as one frame
        frames.append(atoms_rot)

    # write only one trajectory file for this orientation
    write_xyz_trajectory(temp_dir / "trajectory.xyz", frames)

    return temp_dir



def compute_profiles_for_orientation(orient_dir, box_size, bin_size):
    """
    Run volume fraction / number density calculation and compute SLDs.
    """

    # Compute densities from the dvp module
    prev_cwd = os.getcwd()  # baseline directory
    os.chdir(str(orient_dir))  # go to the dir containing xyz

    try:
        dvp.calcAll(size=box_size, w_slice=bin_size)
    finally:
        os.chdir(prev_cwd)

    # Compute SLDs (densities smoothed with sigma=2 Å)
    sld = compute_sld_profiles(orient_dir)

    return sld


def build_lookup_entry(orient_dir, theta, phi, box_size, sld):
    """
    Read output density files for one orientation and build one JSON entry.
    """

    # Build MATLAB-friendly entry
    zC, rhoC, _ = read_density_dat(orient_dir / "density_C.dat")
    _, rhoH, _ = read_density_dat(orient_dir / "density_H.dat")
    _, rhoN, _ = read_density_dat(orient_dir / "density_N.dat")
    _, rhoO, _ = read_density_dat(orient_dir / "density_O.dat")

    nP = []  # add the phosphorus ND if it exists
    dP = orient_dir / "density_P.dat"
    if dP.exists():
        _, nP, _ = read_density_dat(dP)

    zVF, vf, _ = read_density_dat(orient_dir / "vol_frac_along_z.dat")

    # it should correspond though, as in dvp we use the same bin size for volumes and densities
    z = zC if len(zC) == len(vf) else zVF

    entry = {
        "theta": int(theta),
        "phi": int(phi),
        "box_size": box_size,  # added for convenience
        "z": list(map(float, z)),  # convert to a list of floats
        "VF": list(map(float, vf)),
        "nC": list(map(float, rhoC)),
        "nH": list(map(float, rhoH)),
        "nN": list(map(float, rhoN)),
        "nO": list(map(float, rhoO)),
        "nP": list(map(float, nP)) if nP else [],
        "sld_total": list(map(float, sld)),
    }

    return entry


def run_single_orientation(wriggled_pdbs, out_root, theta, phi, box_size, bin_size):
    """
    Run the full calculation for one theta/phi orientation.
    """

    # Create a sub-dir for every combination of theta-phi angles
    orient_dir = out_root / f"theta{theta}_phi{phi}"
    orient_dir.mkdir(parents=True, exist_ok=True)

    # Remove old distributions
    clean_distributions(orient_dir)

    # Rotate structures and write XYZ files
    rotate_wriggled_structures_for_orientation(
        wriggled_pdbs=wriggled_pdbs,
        theta=theta,
        phi=phi,
        orient_dir=orient_dir,
    )

    # Compute densities and SLD
    sld = compute_profiles_for_orientation(
        orient_dir=orient_dir,
        box_size=box_size,
        bin_size=bin_size,
    )

    # Build and return one lookup-table entry
    entry = build_lookup_entry(
        orient_dir=orient_dir,
        theta=theta,
        phi=phi,
        box_size=box_size,
        sld=sld,
    )

    return entry


def build_lookup_table(wriggled_pdbs, out_root, thetas, phis, box_size, bin_size):
    """
    Build the full theta/phi lookup table.

    The full array is a list containing one list for each theta angle.
    Each theta-list contains one dictionary for each phi angle.
    Each dictionary contains something like:

    {
        "theta": theta_angle,
        "phi": phi_angle,
        "VF": [...],
        "nC": [...],
        "nH": [...],
        "nN": [...],
        "nO": [...]
    }
    """

    full_array = []

    for theta in thetas:
        row = []

        for phi in phis:
            entry = run_single_orientation(
                wriggled_pdbs=wriggled_pdbs,
                out_root=out_root,
                theta=theta,
                phi=phi,
                box_size=box_size,
                bin_size=bin_size,
            )

            row.append(entry)

        full_array.append(row)

    return full_array


# ------------------------ main API / pipeline ------------------------

def main_api(
    pdb,
    xyz=None,
    theta_angles="0",
    phi_angles="0",
    outdir="orientations",
    box_size=None,
    bin_size=1.0,
    python=sys.executable,
    skip_wriggle=False,
    wriggle_arg_input="--input",
    wriggle_arg_nstruc="-n",
    nstruc=20,
    json_name="lookup.json",
    wriggle_mode="protein",
    bb_sd=None,
    bb_window=None,
    bb_attempts=None,
    bb_stride=None,
    bb_clash_scale=None,
    sc_sd_backbone=None,
    sc_opt_passes=None,
    sc_opt_scale=None,
    sc_opt_var=None,

):  # this is for MATLAB call with arguments

    root = Path(".").resolve()  # absolute path for output directory
    out_root = root / outdir  # name and location of the output directory
    out_root.mkdir(parents=True, exist_ok=True)  # create the output dir if not exists

    write_run_metadata(
        out_root=out_root,
        pdb=pdb,
        theta_angles=theta_angles,
        phi_angles=phi_angles,
        box_size=box_size,
        bin_size=bin_size,
        nstruc=nstruc,
        skip_wriggle=skip_wriggle,
        json_name=json_name,
        wriggle_mode=wriggle_mode
    )


    # Clean any old distributions in top level
    clean_distributions(root)

    # Angle grids
    thetas = parse_angle_spec(theta_angles)
    phis = parse_angle_spec(phi_angles)

    # Prepare aligned/wriggled structures once
    wriggled_pdbs = prepare_wriggled_structures(
        pdb=pdb,
        out_root=out_root,
        python=python,
        skip_wriggle=skip_wriggle,
        wriggle_arg_input=wriggle_arg_input,
        wriggle_arg_nstruc=wriggle_arg_nstruc,
        nstruc=nstruc,
        wriggle_mode=wriggle_mode,
        bb_sd=bb_sd,
        bb_window=bb_window,
        bb_attempts=bb_attempts,
        bb_stride=bb_stride,
        bb_clash_scale=bb_clash_scale,
        sc_sd_backbone=sc_sd_backbone,
        sc_opt_passes=sc_opt_passes,
        sc_opt_scale=sc_opt_scale,
        sc_opt_var=sc_opt_var,
        )

    # Build lookup table over all theta/phi orientations
    full_array = build_lookup_table(
        wriggled_pdbs=wriggled_pdbs,
        out_root=out_root,
        thetas=thetas,
        phis=phis,
        box_size=box_size,
        bin_size=bin_size,
    )

    # Output list of dicts for MATLAB that can be imported
    json_out = out_root / json_name
    with open(json_out, "w") as final:
        json.dump(full_array, final)

    print("Done. JSON lookup written to", json_out)

    return full_array  # maybe return the results object


# ------------------------ CLI entry point ------------------------

def main(): # this is for CLI call
    arg_parse = argparse.ArgumentParser()
    arg_parse.add_argument("--pdb", default="protein.pdb")
    # New names referring to polar coordinates
    arg_parse.add_argument("--theta-angles", default="0:90:15", help="theta grid in deg, applied about z. Format a:b:step")
    arg_parse.add_argument("--phi-angles", default="0:90:15", help="phi grid in deg, applied about x. Format a:b:step")
    # Old names (aliases)
    arg_parse.add_argument("--x-angles", default="0:90:45", help="DEPRECATED alias for --theta-angles")
    arg_parse.add_argument("--y-angles", default="0:90:45", help="DEPRECATED alias for --phi-angles")
    arg_parse.add_argument("--outdir", default="orientations")
    arg_parse.add_argument("--json", default="lookup.json", help="Output JSON filename (written under --outdir)")
    arg_parse.add_argument("--box-size", type=float, default=60)
    arg_parse.add_argument("--bin", type=float, default=0.5)
    arg_parse.add_argument("--python", default=sys.executable)
    arg_parse.add_argument("--skip-wriggle", action="store_true")
    arg_parse.add_argument("--wriggle-arg-input", default="--input", help="Flag name in wriggle_driver.py that accepts the PDB filename")
    arg_parse.add_argument("--wriggle-arg-nstruc", default="-n", help="Flag name in wriggle_driver.py that accepts the number of structures to generate by wriggling")
    arg_parse.add_argument("--nstruc", default=20, help="Number of wriggled structures")
    arg_parse.add_argument("--wriggle-mode", choices=["protein", "peptide"], default="protein", help="Preset used by wriggle_driver.py for wriggling parameters")

    # Optional wriggling overrides passed through to wriggle_driver.py
    arg_parse.add_argument("--bb-sd", type=float, default=None)
    arg_parse.add_argument("--bb-window", type=int, default=None)
    arg_parse.add_argument("--bb-attempts", type=int, default=None)
    arg_parse.add_argument("--bb-stride", type=int, default=None)
    arg_parse.add_argument("--bb-clash-scale", type=float, default=None)

    arg_parse.add_argument("--sc-sd-backbone", type=float, default=None)
    arg_parse.add_argument("--sc-opt-passes", type=int, default=None)
    arg_parse.add_argument("--sc-opt-scale", type=float, default=None)
    arg_parse.add_argument("--sc-opt-var", type=float, default=None)



    args = arg_parse.parse_args()

    theta_angles = args.theta_angles if args.theta_angles is not None else args.x_angles
    phi_angles = args.phi_angles if args.phi_angles is not None else args.y_angles

    main_api(pdb=args.pdb,  
        theta_angles=theta_angles, phi_angles=phi_angles, outdir=args.outdir,
        box_size=args.box_size, bin_size=args.bin,
        python=args.python, skip_wriggle=args.skip_wriggle, 
        wriggle_mode=args.wriggle_mode,wriggle_arg_input=args.wriggle_arg_input,         wriggle_arg_nstruc=args.wriggle_arg_nstruc,
        nstruc=args.nstruc, json_name=args.json, bb_sd=args.bb_sd,
        bb_window=args.bb_window, bb_attempts=args.bb_attempts,
        bb_stride=args.bb_stride, bb_clash_scale=args.bb_clash_scale,
        sc_sd_backbone=args.sc_sd_backbone, sc_opt_passes=args.sc_opt_passes,
        sc_opt_scale=args.sc_opt_scale, sc_opt_var=args.sc_opt_var,)

if __name__ == "__main__":
    main()

