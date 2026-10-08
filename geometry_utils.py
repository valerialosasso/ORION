import numpy as np

def rotation_matrix_from_vectors(v_from, v_to):
    """
    Function used for initial alignment.
    Return the rotation matrix that rotates vector v_from (= peptide principal axis) onto vector v_to (x-axis = [1, 0, 0]).
    """
    # Normalise both vectors so only direction matters, not magnitude -> dividing by the magnitude
    a = v_from / np.linalg.norm(v_from)
    b = v_to   / np.linalg.norm(v_to)

    # The cross product of the two normalised vectors gives the axis of rotation
    v = np.cross(a, b)

    # The dot product gives cos(theta), where theta is the angle between the two normalised vectors
    c = np.dot(a, b)

    # Length of the cross product = sin(theta)
    s = np.linalg.norm(v)

    # If vectors are already aligned (angle close to 0)
    if s < 1e-8:
        # If dot product is positive → same direction → identity rotation
        if c > 0:
            return np.eye(3) # Return an array with 1s on the diagonal and zeros elsewhere

        # Otherwise they are opposite (180° rotation)
        # We must rotate around any axis perpendicular to "a" and rotate 180 degrees about it to avoid NaNs

        # Pick an arbitrary perpendicular axis in a stable way
        axis = np.array([1, 0, 0])
        if abs(a[0]) > 0.9:
            axis = np.array([0, 1, 0])

        v = np.cross(a, axis)
        v /= np.linalg.norm(v)

        # Rodrigues formula for 180 degrees rotation:
        # R = 2 vv^(T) − I
        return 2 * np.outer(v, v) - np.eye(3)

    # Skew-symmetric cross-product matrix of v
    # This matrix represents "cross with v"
    vx = np.array([
        [ 0,   -v[2],  v[1]],
        [ v[2],  0,   -v[0]],
        [-v[1], v[0],   0 ]
    ])

    # Rodrigues' rotation formula:
    # R = I + vx + vx^2 * ((1 - cosθ) / sin^{2}θ)
    R = np.eye(3) + vx + vx @ vx * ((1 - c) / (s ** 2))

    return R



def center_and_align_atoms_along_x(atoms):
    """
    Center a set of atoms on the origin and align their principal axis along +x.

    Parameters
    ----------
    atoms : (N, 3) ndarray
        Cartesian coordinates of atoms.

    Returns
    -------
    atoms_aligned : (N, 3) ndarray
        Centered and rotated coordinates.

   """

    # ------------------------------------------------------------
    # 1) Center atoms on their centroid
    # ------------------------------------------------------------

    # Compute the geometric center (average position)
    coords = np.array([[x, y, z] for _, x, y, z in atoms], dtype=float)
    center = coords.mean(axis=0)

    # Translate all atoms so the centroid is at (0, 0, 0)
    atoms_centered = coords - center

    # ------------------------------------------------------------
    # 2) Find the principal axis using PCA
    # ------------------------------------------------------------

    # Compute covariance matrix of coordinates
    # This reflects how positions are spread in x, y, z
    cov = np.cov(atoms_centered.T)


    # Decomposition of covariance matrix:
    # Eigenvectors = principal axes
    # Eigenvalues  = variance along those axes
    eigvals, eigvecs = np.linalg.eigh(cov)

    # Eigenvectors are sorted by ascending eigenvalue
    # The largest eigenvalue corresponds to the long axis
    principal_axis = eigvecs[:, np.argmax(eigvals)]

    # ------------------------------------------------------------
    # 3) Rotate principal axis → +x direction
    # ------------------------------------------------------------

    # Define target direction (x-axis)
    x_axis = np.array([1.0, 0.0, 0.0])

    # Compute rotation matrix that aligns principal_axis onto x_axis
    R = rotation_matrix_from_vectors(principal_axis, x_axis)

    # Apply rotation to all atoms

    atoms_aligned = atoms_centered @ R.T

    return atoms_aligned


def center_atoms_on_origin(atoms):
    """Center atoms on origin by subtracting centroid."""
    if not atoms:
        return atoms
    coords = np.array([[x, y, z] for _, x, y, z in atoms], dtype=float)
    centroid = coords.mean(axis=0)
    coords0 = coords - centroid
    return [(atoms[i][0], float(coords0[i, 0]), float(coords0[i, 1]), float(coords0[i, 2]))
            for i in range(len(atoms))]


def rotate_atoms_theta_phi(atoms, theta_deg, phi_deg):
    """
    Rotate a set of atoms using the "protein-axis" definition (spin around it = phi) and its tilt relative to z (theta).

   - theta: a tilt of the protein long axis relative to z.
      This must move atoms in z when theta changes.

    - phi: a spin (twist) around the protein long axis.
      For a long thin protein, spinning around its own long axis should have
         only a small effect on z-profiles.

    This function assumes that before calling it, the
    protein long axis is aligned to the +x direction

    Parameters
    ----------
    atoms : list of tuples like [(line_index, x, y, z), ...]
        - line_index: original PDB line index (so to write coords back)
        - x, y, z: coordinates (already centered)

    theta_deg : float
        - Tilt angle in degrees.

   phi_deg : float
        - Spin angle in degrees.

    Returns
    -------
    rotated_atoms: new list of tuples
        Same format as input: [(line_index, x_rot, y_rot, z_rot), ...]
    """

    # Handle empty input 
    if not atoms:
       return atoms

    # Convert angles from degrees to radians because numpy trigonometric functions expect radians
    th = np.deg2rad(theta_deg)  # theta in radians
    ph = np.deg2rad(phi_deg)    # phi in radians

    # Precompute cos/sin for speed and numerical cleanliness
    cosz, sinz = float(np.cos(th)), float(np.sin(th))   # cos(theta), sin(theta)

    # 1) Apply the "tilt" rotation corresponding to theta 
   # We implement tilt as a rotation about y. This is because:
    #   - the protein long axis is assumed to be aligned with +x.
    #   - Rotating about y swings the +x axis toward +z (or -z), changing z.
    #
    # This is the classic rotation matrix Ry(theta):
    #
    #   [ cos(theta)  0   sin(theta) ]
    #   [  0     1    0   ]
    #   [ -sin(theta)  0   cos(theta) ]
    #
    # After this tilt, the protein long axis will no longer be exactly +x,
    # it will become a tilted direction in the x–z plane.
    Ry = np.array([
        [ cosz, 0.0,  sinz],
        [ 0.0, 1.0,  0.0],
        [-sinz, 0.0,  cosz],
    ], dtype=float)

    # 2) Determine the new protein long axis direction after the tilt
    # Initially, the protein long axis is along +x:
    x_axis = np.array([1.0, 0.0, 0.0], dtype=float)

    # After applying the tilt Ry, that axis becomes:
    axis = Ry @ x_axis
    # axis is a 3-vector giving the long axis direction 

    # Normalise axis to unit length (important for Rodrigues' rotation formula)
    axis_norm = float(np.linalg.norm(axis))
    if axis_norm < 1e-12:
        # protect against divide by zero.
       axis = np.array([1.0, 0.0, 0.0], dtype=float)
    else:
        axis /= axis_norm

    # Unpack axis components
    ux, uy, uz = float(axis[0]), float(axis[1]), float(axis[2])

    # 3) Apply phi as a spin around the protein long axis
    # This is the main difference wrt the earlier implementation
    #
    # Earlier phi meant "rotate about z".
    # But here we want phi to be "rotate about the protein axis".
    #
    # We implement rotation about an arbitrary unit axis using Rodrigues' formula:
    #
    #   R = I cos(phi) + (1-cos(phi)) (u u^T) + sin(phi) [u]_x
    #
    # where:
    #   - I is 3x3 identity matrix
    #   - u is the unit axis vector (ux,uy,uz)
    #   - u u^T is the outer product (a 3x3 matrix)
    #   - [u]_x is the skew-symmetric cross-product matrix:
    #         [  0  -uz   uy ]
    #         [  uz  0   -ux ]
    #         [ -uy  ux   0  ]
    #
    cph, sph = float(np.cos(ph)), float(np.sin(ph))  # cos(phi), sin(phi)

    I = np.eye(3, dtype=float)

    # Skew-symmetric matrix for cross product with axis u
    K = np.array([
        [ 0.0, -uz,  uy],
        [ uz,  0.0, -ux],
        [-uy,  ux,  0.0],
    ], dtype=float)

    # Outer product u u^T
    uuT = np.outer(axis, axis)

    # Rodrigues rotation matrix about axis by phi
    Rphi = I * cph + (1.0 - cph) * uuT + sph * K

    # 4) Combine rotations: tilt first, then spin around new axis
    # We want:
    #   - first apply the theta tilt (Ry)
    #   - then apply the phi spin around the tilted long axis (Rphi)
    #
    # When applying matrices to column vectors, the total is:
    #   v' = Rphi * (Ry * v) = (Rphi * Ry) * v
    #
    # We store the combined matrix as:
    R = Rphi @ Ry

    # 5) Apply rotation to all atom coordinates
    # Extract Nx3 array of coordinates from atoms list
    coords = np.array([[x, y, z] for _, x, y, z in atoms], dtype=float)
   # coords are stored as row vectors in this array.
    # To apply R (which is defined for column vectors), we use coords @ R.T
    coords_r = coords @ R.T

    # Reattach indexes and return the rotated output in the same format as input
    out_atoms = []
    for i in range(len(atoms)):
        line_idx = atoms[i][0]  # original PDB line index (for writing back later)
        out_atoms.append((
            line_idx,
            float(coords_r[i, 0]),
            float(coords_r[i, 1]),
            float(coords_r[i, 2]),
        ))

    return out_atoms


