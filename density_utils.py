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

def compute_sld_profiles(orient_dir):  # read 2 columns from density files
    zC, rhoC, _ = read_density_dat(orient_dir / "density_C.dat")
    zH, rhoH, _ = read_density_dat(orient_dir / "density_H.dat")
    zN, rhoN, _ = read_density_dat(orient_dir / "density_N.dat")
    zO, rhoO, _ = read_density_dat(orient_dir / "density_O.dat")
    z = zC  # values of z

    dz = abs(z[1] - z[0])  # bin size
    sigmaA = 2.0  # hard-coded smoothing sigma in A
    k = gaussian_kernel(dz, sigmaA, truncate=3.0)

    # Smooth number densities
    rhoC_s = convolve_same_length(rhoC, k)
    rhoH_s = convolve_same_length(rhoH, k)
    rhoN_s = convolve_same_length(rhoN, k)
    rhoO_s = convolve_same_length(rhoO, k)

    # SLD
    sldC = [r * B_neutron["C"] for r in rhoC_s]
    sldH = [r * B_neutron["H"] for r in rhoH_s]
    sldN = [r * B_neutron["N"] for r in rhoN_s]
    sldO = [r * B_neutron["O"] for r in rhoO_s]

    # For each z slice, take the SLD contributions from C, H, N, O and add them up into  total SLD
    sldTot = []
    for c, h, n, o in zip(sldC, sldH, sldN, sldO):
        total = c + h + n + o  # sum contributions from all elements
        sldTot.append(total)

    # write output files
    write_two_col(orient_dir / "sld_C.dat", z, sldC)
    write_two_col(orient_dir / "sld_H.dat", z, sldH)
    write_two_col(orient_dir / "sld_N.dat", z, sldN)
    write_two_col(orient_dir / "sld_O.dat", z, sldO)
    write_two_col(orient_dir / "sld_total.dat", z, sldTot)

    # return {
    #   "z": z,              # positions along z-axis (A)
    #  "total": sldTot  # total SLD in A^-2 (not scaled)
    # }
    return sldTot  # total SLD in A^-2 (not scaled)
