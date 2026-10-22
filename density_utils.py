# ------------------------ density I/O helpers ------------------------
import math

#------------------------- constants ---------------------------------

# define all the neutron Bs
B_neutron = {"H": -0.3739e-4, "C":  0.6646e-4, "N":  0.936e-4, "O":  0.5843e-4, "P": 0.513e-4, "D": 0.6671e-4}

def read_density_dat(path):
    zs = []
    means = []
    stds = []
    with open(path) as f:
        for line in f:
            parts = line.split()
            if len(parts) >= 3:
                try:
                    z = float(parts[0])
                    mu = float(parts[1])
                    sd = float(parts[2])
                except ValueError:  # if there is no sd (only 2 columns)
                    continue
                zs.append(z)
                means.append(mu)
                stds.append(sd)
    return zs, means, stds


def write_two_col(path, xs, ys):  # write only z and averages (no sd)
    with open(path, "w") as f:
        for x, y in zip(xs, ys):
            f.write(f"{x:.6f} {y:.12g}\n")


# ------------------------ smoothing and SLD helpers ------------------------



def gaussian_kernel(dz, sigma, truncate=3.0):  # get weights
    half_width_angstrom = truncate * sigma  # physical half-width of kernel in A
    half_width_bins = max(1, int(round(half_width_angstrom / dz)))  # bins covered

    # offsets from -half_width to +half_width in A
    offsets = []
    for bin_index in range(-half_width_bins, half_width_bins + 1):
        offset = bin_index * dz
        offsets.append(offset)

    # non normalised Gaussian weights
    raw_weights = []
    for offset in offsets:
        raw_weight = math.exp(-0.5 * (offset / sigma) ** 2)
        raw_weights.append(raw_weight)

    # normalise so weights sum to 1
    total = sum(raw_weights)
    normalised_weights = [w / total for w in raw_weights]

    return normalised_weights


def convolve_same_length(data, kernel):
    """
    Smooth data by convolving with kernel - .
    Keeps output the same length as input.
    """
    # trivial case: kernel = [1.0] means "no smoothing"
    if len(kernel) == 1:
        return data[:]

    n_data = len(data)
    n_kernel = len(kernel)
    center = n_kernel // 2  # middle index of kernel

    smoothed = []
    for i in range(n_data):
        weighted_sum = 0.0
        for j in range(n_kernel):
            # map kernel index j onto data index
            data_index = i + (j - center)

            # skip kernel weights that fall outside the data range
            if 0 <= data_index < n_data:
                weighted_sum += data[data_index] * kernel[j]

        smoothed.append(weighted_sum)

    return smoothed


# SLD computation (fixed sigma=2 A on densities)

def compute_sld_profiles(orient_dir):
    """Write raw and sigma=2 Å Gaussian-smoothed neutron SLD profiles.

    Include phosphorus when density_P.dat exists. Preserve the historical
    return value: the smoothed total SLD array.
    """
    from pathlib import Path

    orient_dir = Path(orient_dir)
    densities = {}
    z = None
    for element in ("C", "H", "N", "O", "P"):
        path = orient_dir / f"density_{element}.dat"
        if element == "P" and not path.exists():
            continue
        zs, rho, _ = read_density_dat(path)
        if z is None:
            z = zs
        elif len(zs) != len(z) or any(abs(a - b) > 1e-6 for a, b in zip(z, zs)):
            raise ValueError(f"Inconsistent z grid in {path}")
        densities[element] = rho

    if len(z) < 2:
        raise ValueError("At least two z bins are required for SLD calculation")
    dz = abs(z[1] - z[0])
    kernel = gaussian_kernel(dz, sigma=2.0, truncate=3.0)

    raw_total = [sum(densities[e][i] * B_neutron[e] for e in densities)
                 for i in range(len(z))]
    write_two_col(orient_dir / "sld_total_raw.dat", z, raw_total)

    contributions = {}
    for element, rho in densities.items():
        smoothed = convolve_same_length(rho, kernel)
        contributions[element] = [v * B_neutron[element] for v in smoothed]
        write_two_col(orient_dir / f"sld_{element}.dat", z, contributions[element])

    # Avoid retaining an old phosphorus output when the current system has none.
    if "P" not in densities:
        old_p = orient_dir / "sld_P.dat"
        if old_p.exists():
            old_p.unlink()

    smoothed_total = [sum(contributions[e][i] for e in contributions)
                      for i in range(len(z))]
    write_two_col(orient_dir / "sld_total.dat", z, smoothed_total)
    return smoothed_total
