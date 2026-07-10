#!/usr/bin/env python3
"""
stochastic_wriggle.py — peptide/protein torsion "shake" with integrated side‑chain optimiser.

Features
--------
- Side chains: chi rotamer sampling (--sc-sd) + optional quick optimiser (--sc-optimise)
- Backbone: gentle phi/psi tweaks with stride, attempts, and local windows
- Peptide bond omega: small tweaks near trans, with stride/attempts/windows
- Separate clash scales per move type (SC/BB/OMEGA), plus global --clash-scale
- Optional axis realign (to keep or remove overall straightening)
- Reporting of accepted moves

Integrated optimiser 
-------------------------------------------------
  --sc-optimise                enable optimisation by rotamer cleanup after shakes
  --sc-opt-passes N            optimisation passes (default 5)
  --sc-opt-scale X             VDW scale for optimiser (default 0.96; higher = stricter)
  --sc-opt-variation J            variation (deg) around canonical rotamers (default 4.0)
"""

import argparse, math, os, random, sys
from collections import defaultdict, namedtuple

# with namedtuple, every PDB line element can be called by name other than index
Atom = namedtuple("Atom", "atom_index name alternate_location resName chain res_num insertion_code x y z occ b_fact segid element charge")

# ---------- math helpers ----------
def deg2rad(a): 
    return a * math.pi / 180.0 # convert degrees to radians (for angles) 

def rad2deg(a): 
    return a * 180.0 / math.pi # convert radians to degrees (for angles)

def vsub(a,b): # vector subtraction, to compute bond vectors from atom1 to tom2 
    return (a[0]-b[0], a[1]-b[1], a[2]-b[2])

def vdot(a,b): # dot (scalar) product  
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]

def vcross(a,b): # vector perpendicular to both a and b 
    return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])

def vscale(a,s): # scaling a vector 
    return (a[0]*s, a[1]*s, a[2]*s)

def vunit(a): # calculates the unit vector (length 1) 
    n = math.sqrt(max(1e-30, vdot(a,a))) # to avoid division by zero if the vector length is 0
    return (a[0]/n, a[1]/n, a[2]/n), n

def rodrigues_rotate(vec, axis_unit, angle_rad):
    """
    Rotate vector 'vec' around a unit axis 'axis_unit' by 'angle_rad'. Returns the rotated vector. Application of the Rodrigues formula. 
    Every time we rotate,
    1) we compute vec = (atom position − pivot atom position).
    2) We rotate it with this function
    3) We add back the pivot atom position to get the new coordinates of the atom.
    4) We do this for all atoms “downstream” of the bond.
    This way, the whole group of atoms gets twisted.
    """
    ux, uy, uz = axis_unit
    c = math.cos(angle_rad); s = math.sin(angle_rad)
    x, y, z = vec
    return (
        x*(c+ux*ux*(1-c)) + y*(ux*uy*(1-c) - uz*s) + z*(ux*uz*(1-c) + uy*s),
        x*(uy*ux*(1-c) + uz*s) + y*(c+uy*uy*(1-c)) + z*(uy*uz*(1-c) - ux*s),
        x*(uz*ux*(1-c) - uy*s) + y*(uz*uy*(1-c) + ux*s) + z*(c+uz*uz*(1-c)),
    )

def dihedral(p1,p2,p3,p4): # compute torsion angles
    # build bond vectors  
    b0 = vsub(p2,p1) 
    b1 = vsub(p3,p2) # torsion axis
    b2 = vsub(p4,p3)
    b1u, _ = vunit(b1) # normalise torsion axis to unit vector
    v = vsub(b0, vscale(b1u, vdot(b0, b1u))) # projects b0 onto the axis b1u using vdot, then subtracts it out to keep only the perpendicular component. v lies in plane (p1,p2,p3)
    w = vsub(b2, vscale(b1u, vdot(b2, b1u))) # projects b0 onto the axis b1u using vdot,     then subtracts it out to keep only the perpendicular component. w lies in plane (p2,p3,p4).
    # normalise perpendicular components
    v_u, _ = vunit(v)
    w_u, _ = vunit(w)
    x = vdot(v_u, w_u) # cos(angle)
    y = vdot(vcross(b1u, v_u), w_u) # vcross builds a perpendicular vector so we can get the signed angle between the two planes. y = sin(angle) with correct sign 
    return math.atan2(y, x) # angle with a sign

# ---------- PDB I/O ----------
def parse_atom_line(line, ix): # return a namedtuple for each line of the PDB
    rec = line[0:6]
    if rec not in ("ATOM  ", "HETATM"): 
        return None
    try:
        atom_index = int(line[6:11])
        name = line[12:16].strip()
        alternate_location = line[16].strip() 
        resName = line[17:20].strip()
        chain = line[21].strip() or " " 
        res_num = int(line[22:26]) 
        insertion_code = line[26].strip()
        x = float(line[30:38]) 
        y = float(line[38:46]) 
        z = float(line[46:54])
        occ = line[54:60]
        b_fact = line[60:66] 
        segid = line[66:76]
        element = line[76:78].strip()
        charge = line[78:80].strip()
    except Exception:
        return None
    return Atom(atom_index,name,alternate_location,resName,chain,res_num,insertion_code,x,y,z,occ,b_fact,segid,element,charge)

def load_pdb(path):
    """ returns a list of strings (one per line) and a list of Atom namedtuples"""
    with open(path) as f:
        lines = f.readlines()
    atoms = []
    for i, line in enumerate(lines):
        a = parse_atom_line(line, i)
        if a: 
            atoms.append(a)
    return lines, atoms

def write_coords(lines, atoms, out_path):
    """Write PDB output file by matching ATOM/HETATM lines by atom index"""
    coord_map = {a.atom_index: (a.x, a.y, a.z) for a in atoms}

    out = []
    for line in lines:
        rec = line[0:6]
        if rec in ("ATOM  ", "HETATM"):
            try:
                serial = int(line[6:11])
            except Exception:
                out.append(line if line.endswith("\n") else line + "\n")
                continue

            if serial in coord_map:
                x, y, z = coord_map[serial]
                new = list(line.rstrip("\n"))
                if len(new) < 80:
                    new += [" "] * (80 - len(new))
                xyz = f"{x:8.3f}{y:8.3f}{z:8.3f}"
                new[30:54] = list(xyz)
                out.append("".join(new) + "\n")
            else:
                out.append(line if line.endswith("\n") else line + "\n")
        else:
            out.append(line if line.endswith("\n") else line + "\n")

    with open(out_path, "w") as f:
        f.writelines(out)


# ---------- residue / torsion setup ----------

# dictionary that defines the chi angles for the side chains of each amino acid
 
CHI_DEF = {
    "ARG":[("N","CA","CB","CG"),("CA","CB","CG","CD"),("CB","CG","CD","NE"),("CG","CD","NE","CZ")],
    "ASN":[("N","CA","CB","CG"),("CA","CB","CG","OD1")],
    "ASP":[("N","CA","CB","CG"),("CA","CB","CG","OD1")],
    "CYS":[("N","CA","CB","SG")],
    "GLN":[("N","CA","CB","CG"),("CA","CB","CG","CD"),("CB","CG","CD","OE1")],
    "GLU":[("N","CA","CB","CG"),("CA","CB","CG","CD"),("CB","CG","CD","OE1")],
    "HIS":[("N","CA","CB","CG"),("CA","CB","CG","ND1")],
    "ILE":[("N","CA","CB","CG1"),("CA","CB","CG1","CD1")],
    "LEU":[("N","CA","CB","CG"),("CA","CB","CG","CD1")],
    "LYS":[("N","CA","CB","CG"),("CA","CB","CG","CD"),("CB","CG","CD","CE"),("CG","CD","CE","NZ")],
    "MET":[("N","CA","CB","CG"),("CA","CB","CG","SD"),("CB","CG","SD","CE")],
    "PHE":[("N","CA","CB","CG"),("CA","CB","CG","CD1")],
    "PRO":[("N","CA","CB","CG"),("CA","CB","CG","CD")],
    "SER":[("N","CA","CB","OG")],
    "THR":[("N","CA","CB","OG1")],
    "TRP":[("N","CA","CB","CG"),("CA","CB","CG","CD1")],
    "TYR":[("N","CA","CB","CG"),("CA","CB","CG","CD1")],
    "VAL":[("N","CA","CB","CG1")],
}

# set of backbone atoms to check if an atom belongs to backbone
BACKBONE_NAMES = {"N","CA","C","O","OXT","H","HA","H1","H2","H3"}


# ---------- clash model ----------

# dictionary with vdw values
VDW = {"H":1.20, "C":1.70, "N":1.55, "O":1.52, "S":1.80, "P":1.80}

def elem_of(a): # Infer element from atom name's first letter
    nm = a.name.strip() # list of characters of the name
    letters = []
    for ch in nm:
        if ch.isalpha(): # if the character is a letter, append it to letters list
            letters.append(ch)
    base = "".join(letters)
    if base:
       first_letter = base[0].upper()
    else:
       first_letter = "C"
    return first_letter

def vdw_radius(a): # get the vdw radius of the element or take the carbon's as default (1.70) 
    return VDW.get(elem_of(a), 1.70)

def dist2(a,b): # compute squared distance between atoms (for comparisons without sqrt)
    dx=a.x-b.x 
    dy=a.y-b.y 
    dz=a.z-b.z 
    return dx*dx+dy*dy+dz*dz

def same_res(a, b):
    """
    Return True if atoms a and b belong to the same residue.
    We check:
      - chain identifier
      - residue number
      - insertion code
      - residue name
    """
    same_chain = (a.chain == b.chain)
    same_number = (a.res_num == b.res_num)    
    same_icode = (a.insertion_code == b.insertion_code)        # insertion code
    same_resname = (a.resName == b.resName)

    if same_chain and same_number and same_icode and same_resname:
        return True
    else:
        return False

# check if there are clashes between two atoms
def clash_between(a, b, scale=1.0, bond_cut2=1.9*1.9, ignore_same_res=True): # longest single bond = 2.0 for S-S, anything below is safe

    # skip if atoms are in the same residue 
    if ignore_same_res and same_res(a, b):
        return False

    # Squared distance between atom centers
    d2 = dist2(a, b)

    # If they are bonded (distance less than bond cutoff), ignore
    if d2 < bond_cut2:
        return False

    # Effective clash distance = sum of vdw radii * scale (here = 1)
    rsum = (vdw_radius(a) + vdw_radius(b)) * scale

    # If the squared distance is smaller than that threshold → clash
    if d2 < (rsum * rsum):
        return True
    else:
        return False


def any_clash(moved, static, scale, bond_cut2, ignore_same_res):
    """ Loops over every moved atom; compares each to every static atom (all the rest of the structure). If any pair clashes, immediately return True."""
    for a in moved:
        for b in static:
            if clash_between(a,b,scale,bond_cut2,ignore_same_res): 
                return True
    """ Check every atom pair within the moved group """
    for i in range(len(moved)):
        for j in range(i+1,len(moved)):
            if clash_between(moved[i], moved[j], scale, bond_cut2, ignore_same_res): 
                return True
    return False # we reach this point only if no clashes above 

# ---------- transforms ----------

def rotate_group(atoms_subset, pivot_a, pivot_b, angle_deg):
    # define the axis (bond vector from pivot_a → pivot_b, for example CA-CB for chi1)
    axis = vsub((pivot_b.x,pivot_b.y,pivot_b.z),
                (pivot_a.x,pivot_a.y,pivot_a.z))
    axis_u, _ = vunit(axis)   # normalise axis

    # choose the origin at pivot_b, which stays fixed
    origin = (pivot_b.x, pivot_b.y, pivot_b.z)

    # convert angle to radians
    ang = deg2rad(angle_deg)

    # rotate all atoms in atoms_subset
    out = []
    for a in atoms_subset:
        # translate atom so origin is at pivot_b
        v = (a.x - origin[0], a.y - origin[1], a.z - origin[2])
        # rotate vector v around axis by angle ang
        vr = rodrigues_rotate(v, axis_u, ang)
        # translate back to original origin and create new Atom with replaced coordinates
        out.append(a._replace(x=vr[0] + origin[0], # method for (immutable) namedtuples
                              y=vr[1] + origin[1],
                              z=vr[2] + origin[2]))
    return out


def get_chain_axis(all_atoms, chain): # estimate the main axis direction of a chain
    """
    Estimate the axis of a chain by looking at the vector
    from the first Cα atom to the last Cα atom.
    """
    # collect CA atoms for this chain
    cas = []
    for atom in all_atoms:
        if atom.chain == chain and atom.name == "CA":
            cas.append(atom)
    # if fewer than 2, return a default axis
    if len(cas) < 2:
        return (1.0, 0.0, 0.0)
    # sort them by residue number and insertion code
    def sort_key(atom):
        return (atom.res_num, atom.insertion_code)
    cas.sort(key=sort_key)
    # vector from first CA to last CA
    first = cas[0]
    last  = cas[-1]
    v = (last.x - first.x, last.y - first.y, last.z - first.z)
    # normalise to unit vector
    axis, _ = vunit(v)

    return axis


def align_chain_long_axis(all_atoms, chain, u_before_move):
    """
    Rigidly rotate one chain so that its *current* long axis (u_after_move)
    is rotated to match a *target* direction (u_before_move).

    Steps:
      1) Find atoms belonging to this chain.
      2) Compute u_after_move (first to last CA direction).
      3) Rotation axis = u_after_move × u_before_move ; angle = arccos( dot(...) ).
      4) Rotate all chain atoms around the chain COM by that angle.
    """
    # 1) Collect indices of atoms in the target chain
    idxs = []
    for i, a in enumerate(all_atoms):
        if a.chain == chain:
            idxs.append(i)

    # If there aren't enough atoms to define an axis, do nothing
    if len(idxs) < 2:
        return all_atoms

    # 2) Compute the current chain axis (unit vector)
    u_after_move = get_chain_axis(all_atoms, chain)

    # 3) Build rotation axis (cross product of current and target axes)
    axis = vcross(u_after_move, u_before_move)
    axis_u, magnitude = vunit(axis)

    # If the cross product is ~0, the axes are already (anti)parallel
    if magnitude < 1e-8:
        return all_atoms

    # Force the dot to [-1, 1] in case of floating point errors, then angle = arccos(dot)
    dot_val = vdot(u_after_move, u_before_move)
    if dot_val < -1.0:
        dot_val = -1.0
    elif dot_val > 1.0:
        dot_val = 1.0
    ang = math.acos(dot_val)

    # 4) Compute the chain center of mass (simple average of xyz)
    sum_x = 0.0
    sum_y = 0.0
    sum_z = 0.0
    for i in idxs:
        sum_x += all_atoms[i].x
        sum_y += all_atoms[i].y
        sum_z += all_atoms[i].z
    n = float(len(idxs))
    comx = sum_x / n
    comy = sum_y / n
    comz = sum_z / n

    # 5) Rotate each atom of this chain around COM by angle 'ang' about 'axis_u'
    out = list(all_atoms)  # copy
    for i in idxs:
        a = out[i]
        # translate to COM
        vx = a.x - comx
        vy = a.y - comy
        vz = a.z - comz
        # rotate with Rodrigues' formula
        vrx, vry, vrz = rodrigues_rotate((vx, vy, vz), axis_u, ang)
        # translate back and replace coordinates
        out[i] = a._replace(x=vrx + comx, y=vry + comy, z=vrz + comz)

    return out


# ---------- grouping ----------
def group_atoms_by_residue(atoms):
    """
    Group atoms by residue and also build the residue order per chain.

    Returns:
      groups: dict where
              key   = (chain, res_num, insertion_code, resName)
              value = list of Atom objects in that residue (sorted by atom index)
      chains: dict where
              key   = chain ID (e.g., 'A')
              value = list of (res_num, insertion_code, resName) in order along the chain
    """
    # ---- 1) Build groups: (chain, res_num, iCode, resName) -> [atoms]
    groups = {}  # use a plain dict 

    for a in atoms:
        res_num = a.res_num
        i_code  = a.insertion_code
        key = (a.chain, res_num, i_code, a.resName)

        if key not in groups: # not added yet
            groups[key] = [] # create empty list as dict value
        groups[key].append(a) # append atoms entry to list

    # ---- 2) Sort atoms inside each residue by a stable index
    def atom_order_key(atom):
        if hasattr(atom, "atom_index"):
            return atom.atom_index
        else:
            return 0  # last resort

    for key in groups:
        groups[key].sort(key=atom_order_key) # sorting the dict values i.e. the lists of atoms inside the values, not the dictionary itself

    # ---- 3) Build chains: chain -> list of all residue identifiers (with duplicates)
    chains = {}
    for key in groups.keys():
        chain_id, res_num, i_code, res_name = key
        if chain_id not in chains: # not added yet
            chains[chain_id] = [] # add empty list for chains
        chains[chain_id].append((res_num, i_code, res_name)) 

    # ---- 4) For each chain, sort residues and remove duplicates while preserving order
    for chain_id in chains:
        residues = chains[chain_id] # (res_num, i_code, res_name)

        # Sort residues by (res_num, insertion_code); keep resName with them
        pairs = []
        for res in residues:
            res_num, i_code, res_name = res
            sort_key = (res_num, i_code)
            pairs.append((sort_key, res)) # keep resName

        # Sort the pairs by sort_key
        def pair_sort_key(pair):
            # pair is (sort_key, residue)
            # we want to sort by the first element (sort_key)
            return pair[0]

        pairs.sort(key=pair_sort_key)

        # Remove duplicates
        seen = set() # residues we already added 
        ordered = [] # final ordered list without duplicates
        for _, res in pairs: # pairs is (sort_key, res), we ignore the first
            if res not in seen:
                ordered.append(res)
                seen.add(res)

        chains[chain_id] = ordered

    return groups, chains


def atom_lookup(res_atoms):
    """
    Build a dictionary mapping atom names to Atom objects
    for a given residue.

    Example: {"N": Atom(...), "CA": Atom(...), "C": Atom(...), ...}
    """
    lookup = {}
    for atom in res_atoms:
        lookup[atom.name] = atom
    return lookup



# ---------- counters ----------
class MoveCounter:
    """
    Container for counters that track accepted moves
    """
    def __init__(self): 
        self.sc=0 
        self.bb=0 
        self.omega=0
        self.sc_fix=0

# ---------- side chains (chi) sampling ----------
ROTAMERS=[-60.0, 60.0, 180.0] # canonical chi values: gauche-, gauche+, trans

def sidechain_rotamer_pass(
    groups,
    chains,
    sc_sd_deg,
    rng,
    chain_filter,
    clash_mode,
    clash_scale,
    clash_attempts,
    all_atoms_snapshot,
    bond_cut2,
    ignore_same_res,
    mc  # MoveCounter
):
    """
    For every residue with chi torsions:
      - Identify chi torsions from CHI_DEF
      - Build the set of atoms that should also rotate downstream
      - Try to set each chi near canonical rotamers (−60, +60, 180) with small variation.
      - Rotate the "move" atoms and the downstream side-chain atoms about the pivot bond; if no clash, accept first candidate and update coordinates
      - Add the accepted move to counters

    Arguments:
      groups: dict with keys = (chain, res_num, insertion_code, resName) and values = list[Atom] (this residue's atoms)
      chains: dict with keys = chain and values = ordered list of (res_num, insertion_code, resName)
      sc_sd_deg: std dev (degrees) of the Gaussian variation around the chosen rotamer
      rng: random.Random instance
      chain_filter: set of chain IDs to modify (or None for all)
      clash_mode: "off" | "on"
      clash_scale: VDW scale factor for clash detection
      clash_attempts: max tries per chi torsion before giving up
      all_atoms_snapshot: list[Atom] snapshot of the whole structure (to build static atoms)
      bond_cut2: squared bond cutoff 
      ignore_same_res: bool, if True we ignore intra-residue clashes
      mc: MoveCounter object; we increment mc.sc for each accepted χ move

    Returns:
      Updated 'groups' with new coordinates for any accepted moves.
    """

    # Build a lookup: atom_index to Atom for the full structure snapshot.
    # This is used to build the "static" set, i.e. everything not in the "move" group.
    all_by_idx = {}
    for atom in all_atoms_snapshot:
        all_by_idx[atom.atom_index] = atom

    # Iterate across chains
    for chain_id in chains:
        # If user provided a chain filter, skip chains not in it
        if chain_filter is not None and chain_id not in chain_filter:
            continue

        # Iterate residues of this chain in order
        residues_in_chain = chains[chain_id]
        for residue_tuple in residues_in_chain:
            res_num, insertion_code, res_name = residue_tuple
            key = (chain_id, res_num, insertion_code, res_name)

            # Get atoms for this residue
            res_atoms = groups[key]

            # Build a set of atom names present in this residue
            names_in_residue = set()
            for a in res_atoms:
                names_in_residue.add(a.name)

            # Get the chi definitions for this residue type (list of 4-tuples)
            chis = CHI_DEF.get(res_name, [])
            if len(chis) == 0:
                # No chi angles defined for this residue type
                continue

            # Build a fast name to Atom lookup within the residue
            name_to_atom = {}
            for a in res_atoms:
                name_to_atom[a.name] = a

            # Loop over each chi definition (e.g., ("N","CA","CB","CG"))
            for chi_def in chis:
                a_name, b_name, c_name, d_name = chi_def

                # Ensure all four atoms needed for this dihedral exist in this residue
                have_all = (
                    a_name in names_in_residue and
                    b_name in names_in_residue and
                    c_name in names_in_residue and
                    d_name in names_in_residue
                )
                if not have_all:
                    continue

                # Decide which atoms move when rotating around (b - c) axis.
                #  - Exclude backbone atoms
                #  - Exclude the pivot endpoint 'b' (axis point)
                move_atoms = []
                for atom in res_atoms:
                    if atom.name in BACKBONE_NAMES:
                        continue
                    if atom.name == b_name:
                        continue
                    move_atoms.append(atom)

                # Special case for chi1 defined about CA–CB:
                # don't move CB itself (keep the axis endpoint still)
                if (b_name == "CA") and (c_name == "CB"):
                    move_filtered = []
                    for atom in move_atoms:
                        if atom.name == "CB":
                            continue
                        move_filtered.append(atom)
                    move_atoms = move_filtered

                if len(move_atoms) == 0:
                    # Nothing to rotate; skip
                    continue

                # Build the "static" set = all atoms NOT in move_atoms
                move_indexes = set()
                for m in move_atoms:
                    move_indexes.add(m.atom_index)

                static_atoms = []
                for idx, atom_obj in all_by_idx.items():
                    if idx not in move_indexes:
                        static_atoms.append(atom_obj)

                # How many attempts? If clash checks are off, only 1 attempt is needed.
                if clash_mode == "off":
                    attempts = 1
                else:
                    attempts = max(1, clash_attempts)

                accepted = False

                # Repeatedly sample a candidate chi near a canonical rotamer
                # and accept the first non-clashing one.
                for attempt in range(attempts):
                    # Pick a canonical base (−60, +60, 180)
                    base_rotamer = rng.choice(ROTAMERS)

                    # Add a little Gaussian variation (up to sc_sd_deg / 3.0).
                    proposed_angle = base_rotamer + rng.gauss(0.0, sc_sd_deg / 3.0)

                    # Extract the pivot atoms b and c
                    pivot_b = name_to_atom[b_name]
                    pivot_c = name_to_atom[c_name]

                    # Rotate the move group about (b - c) by the *absolute* target angle.
                    rotated_atoms = rotate_group(move_atoms, pivot_b, pivot_c, proposed_angle)
                    # Clash check
                    if clash_mode == "off":
                        clashing = False
                    else:
                        clashing = any_clash(rotated_atoms, static_atoms, clash_scale, bond_cut2, ignore_same_res)

                    if not clashing:
                        # Accept: write updated coordinates back into this residue AND
                        # into the global lookup (so later moves see these coordinates)
                        updated_by_idx = {}
                        for ra in rotated_atoms:
                            updated_by_idx[ra.atom_index] = (ra.x, ra.y, ra.z)

                        new_res_atoms = []
                        for a0 in res_atoms:
                            if a0.atom_index in updated_by_idx:
                                x, y, z = updated_by_idx[a0.atom_index]
                                a1 = a0._replace(x=x, y=y, z=z) # new coordinates
                                new_res_atoms.append(a1)
                                all_by_idx[a1.atom_index] = a1  # keep global snapshot up to date
                            else:
                                new_res_atoms.append(a0) # add it if it doesn't exist

                        res_atoms = new_res_atoms
                        groups[key] = res_atoms
                        mc.sc += 1  # count one accepted side-chain move
                        accepted = True
                        break  # stop trying for this chi; move to the next chi

                # If not accepted, we leave the residue unchanged for this chi and move on.

    # Return updated groups with any accepted chi rotations applied
    return groups






# ---------- backbone phi/psi ----------



def backbone_small_tweaks(
    all_atoms, # flat list of all Atom records (whole structure)
    groups, # dict {(chain, res_num, insertion_code, resName) : [atoms in that residue]}. 
    chains, # dict {chain : ordered list of (res_num, insertion_code, resName)}.
    bb_sd_deg, # standard deviation for phi/psi moves
    rng,
    chain_filter,  
    clash_mode,
    clash_scale,
    bond_cut2,
    ignore_same_res,
    bb_stride, # only modify every K-th residue 
    do_axis_realign, # if True, rigidly rotate the chain afterward so its long axis stays aligned (prevents global drift)
    bb_attempts, # how many times to try a non-clashing phi/psi move before giving up
    bb_window, # how many residues downstream move together when rotating (local window)
    mc  # MoveCounter
):

    """ Makes small, local rotations along the backbone (phi/psi) while checking for clashes 
    psi_i  : rotate downstream of (CA_i - C_i)
    psi_i+1: rotate downstream of (N_{i+1} - CA_{i+1})

    Accept the first non-clashing candidate per torsion up to bb_attempts.

    (Optionally) realign the chain long axis afterward.
    """


    # ---- Build residue-key - atoms lookup for the current coordinates
    # key: (chain, res_num, insertion_code, resName)
    atoms_by_key = {}
    for a in all_atoms:
        key = (a.chain, a.res_num, a.insertion_code, a.resName)
        if key not in atoms_by_key:
            atoms_by_key[key] = []
        atoms_by_key[key].append(a) # value = Atom record

    # We'll progressively update this copy as we accept moves
    out = list(all_atoms)

    # get an atom by name from a given residue key
    def get_atom_from_key(res_key, atom_name):
        if res_key not in atoms_by_key:
            return None
        atoms_here = atoms_by_key[res_key]
        for aa in atoms_here:
            if aa.name == atom_name:
                return aa
        return None

    # Iterate over chains
    for chain_id in chains:
        # Optional chain filter
        if chain_filter is not None and chain_id not in chain_filter:
            continue

        # Record the chain long axis BEFORE moves (for optional realign)
        u_before_move = get_chain_axis(out, chain_id)

        # Build an ordered list of residue keys for this chain
        res_order = []
        residues_in_chain = chains[chain_id]
        for r in residues_in_chain:
            res_num, insertion_code, res_name = r
            res_order.append((chain_id, res_num, insertion_code, res_name))

        # If chain too short, nothing to do
        if len(res_order) < 3:
            continue

        # Walk along residues
        for idx in range(len(res_order)):
            key_i = res_order[idx]

            # Striding: only act on every bb_stride residues
            if bb_stride > 1:
                desired_mod = (bb_stride - 1) % bb_stride
                if (idx % bb_stride) != desired_mod:
                    continue

            # ---------------------------
            # psi_i (CA_i -> C_i)
            # ---------------------------
            CA_i = get_atom_from_key(key_i, "CA")
            C_i  = get_atom_from_key(key_i, "C")

            if (CA_i is not None) and (C_i is not None):
                # Try up to bb_attempts candidates (or 1 if clash checks are off)
                attempts_for_psi = 1 if (clash_mode == "off") else max(1, bb_attempts)

                for attempt in range(attempts_for_psi):
                    # Small random angle (Gaussian around 0)
                    angle_deg = rng.gauss(0.0, bb_sd_deg)
                    angle_rad = deg2rad(angle_deg)

                    # Rotation axis is unit vector along CA_i -> C_i, origin at C_i
                    axis_vec = vsub((C_i.x, C_i.y, C_i.z), (CA_i.x, CA_i.y, CA_i.z))
                    axis_u, _ = vunit(axis_vec)
                    origin = (C_i.x, C_i.y, C_i.z)

                    # Downstream window: residues after i, limited by bb_window
                    if bb_window > 0:
                        downstream_keys = res_order[idx+1 : idx+1+bb_window]
                    else:
                        downstream_keys = res_order[idx+1 : ]

                    # Build the list of atoms to move (all atoms in downstream window)
                    move_atoms = []
                    for rk in downstream_keys:
                        if rk in atoms_by_key:
                            for b in atoms_by_key[rk]:
                                move_atoms.append(b)

                    # If nothing to move, skip psi_i
                    if len(move_atoms) == 0:
                        break  # nothing to do for psi here

                    # Static atoms = everything NOT in move_atoms
                    moving_idx = set()
                    for a in move_atoms:
                        moving_idx.add(a.atom_index)

                    static_atoms = []
                    for a in out:
                        if a.atom_index not in moving_idx:
                            static_atoms.append(a)

                    # Apply Rodrigues' rotation to all move atoms
                    rotated = []
                    for a in move_atoms:
                        v = (a.x - origin[0], a.y - origin[1], a.z - origin[2])
                        vr = rodrigues_rotate(v, axis_u, angle_rad)
                        rotated.append(a._replace(x=vr[0]+origin[0], y=vr[1]+origin[1], z=vr[2]+origin[2])) # replace coordinates 

                    # Clash check
                    if clash_mode != "off": # need to check for clashes
                        if any_clash(rotated, static_atoms, clash_scale, bond_cut2, ignore_same_res):
                            # try another candidate
                            continue

                    # Accept psi_i move: write back into out and atoms_by_key
                    coords_by_idx = {}
                    for a in rotated:
                        coords_by_idx[a.atom_index] = (a.x, a.y, a.z)

                    # Update 'out'
                    new_out = []
                    for a in out:
                        if a.atom_index in coords_by_idx:
                            x, y, z = coords_by_idx[a.atom_index]
                            new_out.append(a._replace(x=x, y=y, z=z))
                        else:
                            new_out.append(a)
                    out = new_out

                    # Update atoms_by_key for all affected downstream residues
                    for rk in downstream_keys:
                        if rk in atoms_by_key:
                            updated_list = []
                            for b in atoms_by_key[rk]:
                                if b.atom_index in coords_by_idx:
                                    x, y, z = coords_by_idx[b.atom_index]
                                    updated_list.append(b._replace(x=x, y=y, z=z))
                                else:
                                    updated_list.append(b)
                            atoms_by_key[rk] = updated_list

                    mc.bb += 1 # increase the counter
                    break  # stop trying ψ_i; proceed to φ_{i+1}

            # ---------------------------
            # phi_{i+1} (N_{i+1} -> CA_{i+1})
            # ---------------------------
            if (idx + 1) < len(res_order): # check that resid i + 1 exists
                # get the residue-key tuple for residue i+1:
                # key format is (chain_id, res_num, insertion_code, res_name)
                key_ip1 = res_order[idx + 1]
                # get N and CA (forming phi) from resid i + 1 
                N_ip1  = get_atom_from_key(key_ip1, "N")
                CA_ip1 = get_atom_from_key(key_ip1, "CA")

                if (N_ip1 is not None) and (CA_ip1 is not None): # if both atom exist
                    attempts_for_phi = 1 if (clash_mode == "off") else max(1, bb_attempts)

                    for attempt in range(attempts_for_phi):
                        angle_deg = rng.gauss(0.0, bb_sd_deg) # pick a small random angle
                        angle_rad = deg2rad(angle_deg)

                        # set the rotation axis = unit vector from N_{i+1} to CA_{i+1}
                        # origin at CA_{i+1}
                        axis_vec = vsub((CA_ip1.x, CA_ip1.y, CA_ip1.z), (N_ip1.x, N_ip1.y, N_ip1.z))
                        axis_u, _ = vunit(axis_vec)
                        origin = (CA_ip1.x, CA_ip1.y, CA_ip1.z)

                        # rotate a downstream window from residue i+1, around the rotation axis
                        if bb_window > 0:
                            downstream_keys = res_order[idx+1 : idx+1+bb_window]
                        else:
                            downstream_keys = res_order[idx+1 : ]

                        # Move atoms = all atoms in downstream window EXCEPT the N atom itself
                        move_atoms = []
                        for rk in downstream_keys:
                            if rk in atoms_by_key:
                                for b in atoms_by_key[rk]:
                                    if b.name != "N":
                                        move_atoms.append(b)

                        if len(move_atoms) == 0:
                            break  # nothing to do for phi here

                        # Static atoms
                        moving_idx = set()
                        for a in move_atoms:
                            moving_idx.add(a.atom_index)

                        static_atoms = []
                        for a in out:
                            if a.atom_index not in moving_idx:
                                static_atoms.append(a)

                        # Rotate
                        rotated = []
                        for a in move_atoms:
                            v = (a.x - origin[0], a.y - origin[1], a.z - origin[2])
                            vr = rodrigues_rotate(v, axis_u, angle_rad)
                            rotated.append(a._replace(x=vr[0]+origin[0], y=vr[1]+origin[1], z=vr[2]+origin[2]))

                        # Check clashes
                        if clash_mode != "off":
                            if any_clash(rotated, static_atoms, clash_scale, bond_cut2, ignore_same_res):
                                continue

                        # Accept phi_i+1 move
                        coords_by_idx = {}
                        for a in rotated:
                            coords_by_idx[a.atom_index] = (a.x, a.y, a.z)

                        new_out = []
                        for a in out:
                            if a.atom_index in coords_by_idx:
                                x, y, z = coords_by_idx[a.atom_index]
                                new_out.append(a._replace(x=x, y=y, z=z))
                            else:
                                new_out.append(a)
                        out = new_out
                        # update downstream residues
                        for rk in downstream_keys:
                            if rk in atoms_by_key:
                                updated_list = []
                                for b in atoms_by_key[rk]:
                                    if b.atom_index in coords_by_idx:
                                        x, y, z = coords_by_idx[b.atom_index]
                                        updated_list.append(b._replace(x=x, y=y, z=z))
                                    else:
                                        updated_list.append(b)
                                atoms_by_key[rk] = updated_list

                        mc.bb += 1 # update counter
                        break  # done with φ_{i+1}

        # Optional: keep the overall chain pointing the same way as before
        if do_axis_realign:
            out = align_chain_long_axis(out, chain_id, u_before_move)

    return out




# ---------- omega tweaks ----------
def backbone_omega_tweaks(
    all_atoms,
    groups,
    chains,
    om_sd_deg,         # std dev (deg) for small omega variations
    om_max_dev_deg,    # max allowed deviation from trans (in degrees)
    rng,               # random.Random
    chain_filter,      # None or set of chain IDs to include
    clash_mode,        # "off" | "on"
    clash_scale,       # VDW scale for clash check (smaller = more permissive)
    bond_cut2,         # squared bond cutoff (Å^2) for clash_between
    ignore_same_res,   # usually True for speed
    omega_stride,      # act on every k-th peptide (1 = every one)
    no_pro_omega,      # if True, skip peptide bonds leading into PRO residues
    omega_attempts,    # retries per peptide
    omega_window,      # how many residues downstream to move (local bend)
    mc                 # MoveCounter, we increment mc.omega on accept
):
    """
    Small peptide-bond (omega) tweaks:
      - Bond axis: C_i - N_{i+1}
      - Rotate a small downstream window of residues by a few degrees (Gaussian)
      - Reject if clashes or if omega drifts too far from trans (~180) or near cis (~0)
    """

    # --- Build a residue-key -> atoms lookup for current positions
    # Key format: (chain, res_num, insertion_code, resName)
    atoms_by_key = {}
    for a in all_atoms:
        key = (a.chain, a.res_num, a.insertion_code, a.resName)
        if key not in atoms_by_key:
            atoms_by_key[key] = []
        atoms_by_key[key].append(a)

    # We'll progressively update this copy as we accept moves
    out = list(all_atoms)

    # get an atom by name from a given residue key
    def get_atom(res_key, atom_name):
        if res_key not in atoms_by_key:
            return None
        for a in atoms_by_key[res_key]:
            if a.name == atom_name:
                return a
        return None

    # ---- Iterate chains
    for chain_id in chains:
        # Optional chain filter
        if chain_filter is not None and chain_id not in chain_filter:
            continue

        # Ordered residue keys for this chain
        res_order = []
        for (res_num, insertion_code, res_name) in chains[chain_id]:
            res_order.append((chain_id, res_num, insertion_code, res_name))

        # Walk through peptide bonds: between res[i-1] and res[i]
        for idx in range(1, len(res_order)):
            # Stride: only act on every 'omega_stride'-th peptide
            if omega_stride > 1:
                if (idx % omega_stride) != 0:
                    continue

            key_prev = res_order[idx - 1]   # residue i-1 (for C_i-1)
            key_curr = res_order[idx]       # residue i (for N_i, CA_i, C_i)

            # Optionally skip if current residue is PRO
            # (omega in prolines is constrained)
            if no_pro_omega and key_curr[3] == "PRO":
                continue

            # Get the four atoms needed:
            #   C_prev = C of residue i-1
            #   N_curr, CA_curr, C_curr = N/CA/C of residue i
            C_prev = get_atom(key_prev, "C")
            N_curr  = get_atom(key_curr, "N")
            CA_curr = get_atom(key_curr, "CA")
            C_curr  = get_atom(key_curr, "C")

            # If any are missing, skip safely
            if (C_prev is None) or (N_curr is None) or (CA_curr is None) or (C_curr is None):
                continue

            # How many tries for this peptide?
            attempts = 1 if (clash_mode == "off") else max(1, omega_attempts)

            # Try several small variations; accept the first valid one
            accepted = False
            for attempt in range(attempts):
                # Small random delta (deg) for omega
                delta_deg = rng.gauss(0.0, om_sd_deg)
                delta_rad = deg2rad(delta_deg)

                # Rotation axis for omega is along the peptide bond C_{i} - N_{i+1}.
                # Here we use C_prev (of res[i-1]) and N_curr (of res[i]) because we are
                # iterating idx as "current residue"
                axis_vec = (N_curr.x - C_prev.x, N_curr.y - C_prev.y, N_curr.z - C_prev.z)
                axis_u, _ = vunit(axis_vec)

                # Rotate around N_curr (anchor at N_curr)
                origin = (N_curr.x, N_curr.y, N_curr.z)

                # Downstream window to move: current residue onward
                if omega_window > 0:
                    downstream_keys = res_order[idx : idx + omega_window]
                else:
                    downstream_keys = res_order[idx : ]

                # Build list of atoms to move (all atoms in downstream window)
                move_atoms = []
                for rk in downstream_keys:
                    if rk in atoms_by_key:
                        for b in atoms_by_key[rk]:
                            move_atoms.append(b)

                if len(move_atoms) == 0:
                    break  # nothing to move here

                # Static atoms = everything *not* in move_atoms
                moving_idx = set()
                for a in move_atoms:
                    moving_idx.add(a.atom_index)

                static_atoms = []
                for a in out:
                    if a.atom_index not in moving_idx:
                        static_atoms.append(a)

                # Apply Rodrigues' rotation to all move atoms by delta
                rotated = []
                for a in move_atoms:
                    v = (a.x - origin[0], a.y - origin[1], a.z - origin[2])
                    vr = rodrigues_rotate(v, axis_u, delta_rad)
                    rotated.append(a._replace(x=vr[0] + origin[0],
                                              y=vr[1] + origin[1],
                                              z=vr[2] + origin[2]))

                # Clash gate
                if clash_mode != "off":
                    if any_clash(rotated, static_atoms, clash_scale, bond_cut2, ignore_same_res):
                        # Try another delta
                        continue

                # ---- Omega gate: keep near trans, avoid cis
                # Map the potentially moved atoms back by index so we can read
                # their updated coordinates for omega computation.
                rotated_by_idx = {}
                # build a dictionary that maps atom_index → Atom object (with updated coords)
                for a in rotated:
                    rotated_by_idx[a.atom_index] = a

                # Try to get the rotated N atom; if not rotated, keep the original
                if N_curr.atom_index in rotated_by_idx:
                    Nn = rotated_by_idx[N_curr.atom_index]
                else:
                    Nn = N_curr

                # Same for CA
                if CA_curr.atom_index in rotated_by_idx:
                    CAn = rotated_by_idx[CA_curr.atom_index]
                else:
                    CAn = CA_curr

                # Same for C
                if C_curr.atom_index in rotated_by_idx:
                    Cn = rotated_by_idx[C_curr.atom_index]
                else:
                    Cn = C_curr

                # Compute new omega = dihedral(C_prev, N_curr, CA_curr, C_curr)
                omega_rad = dihedral(
                    (C_prev.x, C_prev.y, C_prev.z),
                    (Nn.x,   Nn.y,   Nn.z),
                    (CAn.x,  CAn.y,  CAn.z),
                    (Cn.x,   Cn.y,   Cn.z),
                )
                omega_deg = rad2deg(omega_rad)

                # Normalise dev from trans (180): bring omega into the range (0..360), then
                # compute distance to 180; require dev <= om_max_dev_deg
                omega_mod = (omega_deg + 360.0) % 360.0
                dev_from_trans = abs(omega_mod - 180.0)

                if (dev_from_trans > om_max_dev_deg) or (dev_from_trans < 0.0):
                    # too far from trans — reject
                    continue

                # Avoid cis region (|ω| ~ 0). 
                if (abs(omega_deg) < 150.0) or (abs(omega_deg - 360.0) < 210.0):
                    # too close to cis-like values — reject
                    continue

                # ---- Accept: write back coordinates into 'out' and 'atoms_by_key'
                coords_by_idx = {}
                for a in rotated:
                    coords_by_idx[a.atom_index] = (a.x, a.y, a.z)

                # Update flat list 'out'
                new_out = []
                for a in out:
                    if a.atom_index in coords_by_idx:
                        x, y, z = coords_by_idx[a.atom_index]
                        new_out.append(a._replace(x=x, y=y, z=z))
                    else:
                        new_out.append(a)
                out = new_out

                # Update residue-wise storage 'atoms_by_key'
                for rk in downstream_keys:
                    if rk in atoms_by_key:
                        updated_list = []
                        for b in atoms_by_key[rk]:
                            if b.atom_index in coords_by_idx:
                                x, y, z = coords_by_idx[b.atom_index]
                                updated_list.append(b._replace(x=x, y=y, z=z))
                            else:
                                updated_list.append(b)
                        atoms_by_key[rk] = updated_list

                mc.omega += 1
                accepted = True
                break  # done with this peptide bond

            # If not accepted after all attempts, leave coordinates unchanged for this peptide

    return out


# ---------- side‑chain optimiser ----------

def overlap_score_pair(a, b, scale, bond_cut2, ignore_same_res=True):
    """
    Return a non-negative penalty score for how much two atoms overlap.
    0.0 means no penalty (no clash).
    """
    # Optionally ignore pairs inside the same residue
    if ignore_same_res and same_res(a, b):
        return 0.0

    # Squared distance between atoms
    d2 = dist2(a, b)

    # If within bond cutoff, assume bonded (or very near); do not penalize
    if d2 < bond_cut2:
        return 0.0

    # Effective "nonbonded contact" distance threshold
    rsum = (vdw_radius(a) + vdw_radius(b)) * scale
    rsum2 = rsum * rsum

    # If outside or at VDW sum → no overlap → no penalty
    if d2 >= rsum2:
        return 0.0

    # Otherwise, they overlap: penalty grows with the depth of penetration
    # Use sqrt(d2) but protect against tiny negatives from float noise
    d = math.sqrt(max(1e-12, d2))
    overlap = rsum - d
    return overlap * overlap

def residue_overlap_score(moved_atoms, static_atoms, scale, bond_cut2):
    """
    Sum of pair penalties:
      - moved vs static (groups)
      - moved internal pairs (within the moved group)
    """
    score = 0.0

    # moved vs static
    for a in moved_atoms:
        for b in static_atoms:
            score += overlap_score_pair(a, b, scale, bond_cut2, ignore_same_res=True)

    # moved internal pairs (unique pairs only)
    for i in range(len(moved_atoms)):
        ai = moved_atoms[i]
        for j in range(i + 1, len(moved_atoms)):
            aj = moved_atoms[j]
            score += overlap_score_pair(ai, aj, scale, bond_cut2, ignore_same_res=True)

    return score

def sc_quick_optimiser(
    atoms,          # flat list of all Atom objects
    passes,         # how many passes (stop earlier if no changes)
    scale,          # VDW scale for scoring (e.g., 0.97–0.99 stricter; 0.92–0.95 looser)
    variation,      # Gaussian variation (deg) around canonical rotamers; 0.0 = exact rotamers
    bond_cut,       # bond cutoff (Å)
    seed,           # RNG seed for reproducibility
    chains_filter,  # None or set of chain IDs to include
    report,         # bool; if True, print per-pass change count to stderr
    chi_def=CHI_DEF # dictionary of chi definitions per residue type
):
    rng = random.Random(seed)

    # Group atoms by residue and also get per-chain residue order
    groups, chains = group_atoms_by_residue(atoms)  # groups[(chain,res_num,iCode,resName)] -> [Atom]
    keys = list(groups.keys())                      # keep a stable list of residue keys

    bond_cut2 = bond_cut * bond_cut
    total_changes = 0

    # Passes: keep improving until no more gains (or pass limit reached)
    for p in range(passes):
        changes_this_pass = 0

        for key in keys:
            chain_id, res_num, i_code, res_name = key

            # Optional chain filter
            if (chains_filter is not None) and (chain_id not in chains_filter):
                continue

            res_atoms = groups[key]

            # Build set of names present and a name->Atom map
            names_in_res = set()
            for a in res_atoms:
                names_in_res.add(a.name)

            lu = atom_lookup(res_atoms)   # fast lookup by atom name within residue

            # Skip residues with no chi defined
            chis = chi_def.get(res_name, [])
            if len(chis) == 0:
                continue

            improved_here = False

            # Build static atom list = all atoms not in this residue
            static_atoms = []
            for k2, atoms2 in groups.items():
                if k2 == key:
                    continue
                for a2 in atoms2:
                    static_atoms.append(a2)

            # Consider each chi torsion for this residue
            for chi in chis:
                a_nm, b_nm, c_nm, d_nm = chi

                # Ensure requisite atoms exist
                if not (a_nm in names_in_res and b_nm in names_in_res and
                        c_nm in names_in_res and d_nm in names_in_res):
                    continue

                # Decide which atoms move when rotating around (b -> c)
                move_atoms = []
                for x in res_atoms:
                    if x.name in BACKBONE_NAMES:
                        continue
                    if x.name == b_nm:
                        continue
                    move_atoms.append(x)

                # Special case for chi1 axis CA->CB: do not move CB itself
                if (b_nm == "CA") and (c_nm == "CB"):
                    filtered = []
                    for x in move_atoms:
                        if x.name != "CB":
                            filtered.append(x)
                    move_atoms = filtered

                if len(move_atoms) == 0:
                    continue

                # Compute current chi angle (deg)
                cur_rad = dihedral(
                    (lu[a_nm].x, lu[a_nm].y, lu[a_nm].z),
                    (lu[b_nm].x, lu[b_nm].y, lu[b_nm].z),
                    (lu[c_nm].x, lu[c_nm].y, lu[c_nm].z),
                    (lu[d_nm].x, lu[d_nm].y, lu[d_nm].z)
                )
                cur_deg = rad2deg(cur_rad)

                # Score of the current conformation
                current_score = residue_overlap_score(move_atoms, static_atoms, scale, bond_cut2)

                # Build candidate list (score, target_theta, rotated_atoms)
                candidates = []

                # Try each canonical rotamer ± variation
                for base in ROTAMERS:  # [-60, +60, 180]
                    if variation > 0.0:
                        target_theta = base + rng.gauss(0.0, variation)
                    else:
                        target_theta = base

                    # Rotate by the DELTA needed to reach target from current
                    delta_deg = target_theta - cur_deg

                    rotated = rotate_group(move_atoms, lu[b_nm], lu[c_nm], delta_deg)

                    sc = residue_overlap_score(rotated, static_atoms, scale, bond_cut2)

                    # Store candidate
                    candidates.append((sc, target_theta, rotated))

                # Pick the candidate with the lowest score
                best_sc = None
                best_theta = None
                best_rotated = None
                for tuple_rot in candidates:
                    sc, theta, rot = tuple_rot
                    if (best_sc is None) or (sc < best_sc):
                        best_sc = sc
                        best_theta = theta
                        best_rotated = rot

                # If better than current, accept the best candidate
                # (tiny epsilon to avoid floating error)
                if best_sc + 1e-9 < current_score:
                    # Update residue atoms with new coordinates
                    moved_map = {}
                    for a in best_rotated:
                        moved_map[a.atom_index] = (a.x, a.y, a.z)

                    new_res_atoms = []
                    for a0 in res_atoms:
                        if a0.atom_index in moved_map:
                            x, y, z = moved_map[a0.atom_index]
                            new_res_atoms.append(a0._replace(x=x, y=y, z=z))
                        else:
                            new_res_atoms.append(a0)

                    res_atoms = new_res_atoms
                    groups[key] = res_atoms
                    lu = atom_lookup(res_atoms)  # update lookup after modification

                    improved_here = True
                    changes_this_pass += 1

            # Save any improvements for this residue
            if improved_here:
                groups[key] = res_atoms

        total_changes += changes_this_pass

        if report:
            sys.stderr.write(f"[sc-opt pass {p+1}] changes: {changes_this_pass}\n")

        # Early stop if no changes this pass
        if changes_this_pass == 0:
            break

    # Reassemble flat atom list in original line order
    by_line = {}
    for a in atoms:
        by_line[a.atom_index] = a

    for _, res_atoms in groups.items():
        for a in res_atoms:
            by_line[a.atom_index] = a

    atoms_out = []
    for i in sorted(by_line.keys()):
        atoms_out.append(by_line[i])

    return atoms_out, total_changes



# ---------- MAIN ----------
def main():
    arg_parser=argparse.ArgumentParser(description="Torsion-aware shake with chi/phi/psi and omega tweaks + side‑chain optimiser")
    arg_parser.add_argument("inp")
    arg_parser.add_argument("-o","--out",default=None)
    arg_parser.add_argument("--sc-sd",type=float,default=20.0)
    arg_parser.add_argument("--bb-sd",type=float,default=2.0)
    arg_parser.add_argument("--bb-stride",type=int,default=3, help="apply phi/psi every K residues (1 = every residue)")
    arg_parser.add_argument("--omega-sd",type=float,default=1.0)
    arg_parser.add_argument("--omega-max-dev",type=float,default=10.0)
    arg_parser.add_argument("--omega-stride",type=int,default=1)
    arg_parser.add_argument("--no-pro-omega",default=False)
    arg_parser.add_argument("--seed",type=int,default=None)
    arg_parser.add_argument("--chains",type=str,default=None)
    arg_parser.add_argument("--no-backbone",action="store_true")
    arg_parser.add_argument("--no-omega",action="store_true")
    arg_parser.add_argument("--no-axis-realign",action="store_true", help="do not re-align chain long axis")
    arg_parser.add_argument("--no-sc", action="store_true", help="skip initial sidechain rotamer shake stage")

    # clash controls
    arg_parser.add_argument("--clash-mode",choices=["off","on"],default="on")
    arg_parser.add_argument("--clash-scale",type=float,default=0.85)
    arg_parser.add_argument("--sc-clash-scale",type=float,default=None)
    arg_parser.add_argument("--bb-clash-scale",type=float,default=None)
    arg_parser.add_argument("--omega-clash-scale",type=float,default=None)
    arg_parser.add_argument("--clash-attempts",type=int,default=10)
    arg_parser.add_argument("--bond-cut",type=float,default=1.9)
    arg_parser.add_argument("--no-same-res-clash",default=False)

    # attempts/windows
    arg_parser.add_argument("--bb-attempts",type=int,default=20)
    arg_parser.add_argument("--bb-window",type=int,default=4, help="downstream residues moved by a backbone tweak (0 = full)")
    arg_parser.add_argument("--omega-attempts",type=int,default=20)
    arg_parser.add_argument("--omega-window",type=int,default=4, help="downstream residues moved by an omega tweak (0 = full)")

    # integrated side‑chain optimiser
    arg_parser.add_argument("--sc-optimise",action="store_true", help="run quick side-chain optimiser after shakes")
    arg_parser.add_argument("--sc-opt-passes",type=int,default=5)
    arg_parser.add_argument("--sc-opt-scale",type=float,default=0.96)
    arg_parser.add_argument("--sc-opt-variation",type=float,default=4.0)

    arg_parser.add_argument("--report",action="store_true")
    args=arg_parser.parse_args()

    # seed, input, output
    rng=random.Random(args.seed) if args.seed is not None else random.Random()
    if not os.path.isfile(args.inp):
        sys.stderr.write(f"Error: input PDB not found: {args.inp}\n"); sys.exit(1)
    out_path=args.out or os.path.splitext(args.inp)[0]+"_torsion_shaken.pdb"

    # process PDB
    lines,atoms=load_pdb(args.inp)
    groups,chains=group_atoms_by_residue(atoms)
    chain_filter=set(args.chains.split(",")) if args.chains else None
    bond_cut2=args.bond_cut*args.bond_cut
    ignore_same_res = not args.no_same_res_clash

    # resolve clash scales
    sc_scale   = args.sc_clash_scale   if args.sc_clash_scale   is not None else args.clash_scale
    bb_scale   = args.bb_clash_scale   if args.bb_clash_scale   is not None else args.clash_scale
    omega_scale= args.omega_clash_scale if args.omega_clash_scale is not None else args.clash_scale
    sc_opt_scale = args.sc_opt_scale

    mc=MoveCounter()

    """" Start process """

    # 1. Side chains (initial sampling)
    # Use --no-sc to skip the rotamer "shake" stage (recommended for globular proteins).
    if args.no_sc:
        atoms_sc = atoms
    else:
        groups=sidechain_rotamer_pass(groups,chains,args.sc_sd,rng,chain_filter,
                                      args.clash_mode,sc_scale,args.clash_attempts,
                                      atoms,bond_cut2,ignore_same_res, mc)
        # Flatten residue groups back into a single atom list, in original order
        # dictionary mapping atom_index -> original atom
        by_line = {a.atom_index: a for a in atoms}
        for _, res_atoms in groups.items():
            for a in res_atoms:
                by_line[a.atom_index] = a
        # rebuild a flat list, sorted by line number (restores PDB order)
        atoms_sc = []
        for line_num in sorted(by_line.keys()):
            atoms_sc.append(by_line[line_num])

    # 2. Backbone phi/psi
    if args.no_backbone:
        atoms_bb=atoms_sc
    else:
        atoms_bb=backbone_small_tweaks(atoms_sc,groups,chains,args.bb_sd,rng,chain_filter,
                                       args.clash_mode,bb_scale,bond_cut2,ignore_same_res,
                                       max(1,args.bb_stride), not args.no_axis_realign,
                                       max(1,args.bb_attempts), max(0,args.bb_window), mc)

    # 3. Backbone omega
    if args.no_omega:
        atoms_new=atoms_bb
    else:
        atoms_new=backbone_omega_tweaks(atoms_bb,groups,chains,args.omega_sd,args.omega_max_dev,rng,chain_filter,
                                        args.clash_mode,omega_scale,bond_cut2,ignore_same_res,
                                        max(1,args.omega_stride), args.no_pro_omega,
                                        max(1,args.omega_attempts), max(0,args.omega_window), mc)

    # 4. Side‑chain optimiser 
    if args.sc_optimise:
        # rebuild groups from current atoms so we don't overwrite backbone/omega
        groups, _ = group_atoms_by_residue(atoms_new) # we don't need SC info
        atoms_opt, nchanges = sc_quick_optimiser( # optimiser also returns tot of changes~)
            atoms_new,
            passes=max(1,args.sc_opt_passes),
            scale=args.sc_opt_scale,
            variation=max(0.0,args.sc_opt_variation),
            bond_cut=args.bond_cut,
            seed=args.seed if args.seed is not None else 1,
            chains_filter=chain_filter,
            report=args.report
        )
        atoms_new = atoms_opt
        if args.report:
            sys.stderr.write(f"[report] sc-opt total changes: {nchanges}\n")

    # 5. writing and report
    write_coords(lines,atoms_new,out_path)
    if args.report:
        sys.stderr.write(f"[report] accepted: SC={mc.sc}, BB={mc.bb}, OMEGA={mc.omega}\n")
    sys.stderr.write(f"Done. Wrote {out_path}\n")

if __name__=="__main__":
    main()
