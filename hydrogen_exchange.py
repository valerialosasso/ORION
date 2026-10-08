"""Topology-file-free H classification for ORION's solvent-exposure model.

Coordinates and distances are in Angstrom; SASA is returned in Angstrom squared.
This is a binary accessibility model, not a prediction of exchange kinetics.
"""
import numpy as np
from scipy.spatial import cKDTree

# Conservative covalent H-parent distance windows, not H-bond cutoffs.
PARENT_MAX_A = {'C': 1.25, 'N': 1.20, 'O': 1.15, 'S': 1.50}


def infer_hydrogen_parents(atoms):
    """Return H index -> unique C/N/O/S parent, or fail on ambiguous geometry.

    No atom-name convention or force-field topology is required. The input must
    be a complete, non-periodically-split structure with explicit hydrogens.
    """
    elements = np.array([a[0] for a in atoms])
    xyz = np.array([a[1:4] for a in atoms], dtype=float)
    if not np.isfinite(xyz).all():
        raise ValueError('Non-finite atomic coordinates.')
    hs = np.flatnonzero(elements == 'H')
    if not len(hs):
        raise ValueError('Hydrogen splitting requires a PDB with explicit hydrogens; none found.')
    heavy = np.flatnonzero(elements != 'H')
    tree = cKDTree(xyz[heavy])
    parents = {}
    for h in hs:
        candidates = []
        for local in tree.query_ball_point(xyz[h], max(PARENT_MAX_A.values())):
            p = heavy[local]
            distance = np.linalg.norm(xyz[h] - xyz[p])
            if 0.60 <= distance <= PARENT_MAX_A.get(elements[p], 0):
                candidates.append(int(p))
        if len(candidates) != 1:
            raise ValueError(
                f'H atom {h + 1}: expected one covalent C/N/O/S parent, found '
                f'{len(candidates)} ({[p + 1 for p in candidates]}). '
                'Check protonation, clashes, alternate locations and coordinate integrity.'
            )
        parents[int(h)] = candidates[0]
    return parents


def classify_hydrogens(atoms, sasa_cutoff_A2=0.0, probe_radius_A=1.4,
                       n_sphere_points=960):
    """Return per-atom group-2 mask and per-H audit records.

    All supplied atoms occlude solvent. Supply solute only (no explicit water).
    Group 1 is all remaining H. Heavy atoms are never members of either group.
    Zero cutoff means any sampled accessible surface; it is not a calibrated
    experimental exchange threshold. MDTraj uses elemental van der Waals radii.
    """
    if not np.isfinite(sasa_cutoff_A2) or sasa_cutoff_A2 < 0:
        raise ValueError('SASA cutoff must be finite and nonnegative (A^2).')
    if not np.isfinite(probe_radius_A) or probe_radius_A <= 0:
        raise ValueError('Probe radius must be finite and positive (A).')
    if int(n_sphere_points) != n_sphere_points or n_sphere_points < 1:
        raise ValueError('n_sphere_points must be a positive integer.')
    parents = infer_hydrogen_parents(atoms)
    try:
        import mdtraj as md
    except ImportError as exc:
        raise ImportError('Hydrogen splitting needs MDTraj: pip install mdtraj') from exc
    # SASA only needs elements and coordinates, not a bonded topology file.
    topology = md.Topology()
    residue = topology.add_residue('SOL', topology.add_chain())
    for i, atom in enumerate(atoms):
        topology.add_atom(f'{atom[0]}{i}', md.element.get_by_symbol(atom[0]), residue)
    xyz_nm = np.array([a[1:4] for a in atoms], dtype=float)[None, :, :] / 10.0
    trajectory = md.Trajectory(xyz_nm, topology)
    sasa = md.shrake_rupley(trajectory, mode='atom',
                           probe_radius=probe_radius_A / 10.0,
                           n_sphere_points=int(n_sphere_points))[0] * 100.0
    group2 = np.zeros(len(atoms), dtype=bool)
    records = []
    for h, p in parents.items():
        labile = atoms[p][0] in {'N', 'O', 'S'}
        group2[h] = labile and sasa[h] > sasa_cutoff_A2
        records.append({'atom_index_1based': h + 1, 'parent_index_1based': p + 1,
                        'parent_element': atoms[p][0], 'sasa_A2': float(sasa[h]),
                        'group': 2 if group2[h] else 1})
    return group2, records


def classify_pdbs(paths, **options):
    """Classify each conformation before rotation; reject water/multiple models.

    Atom order must be preserved across the ensemble, as in ORION wriggling.
    PDB names/serials are retained in the audit report, not needed for bonding.
    """
    from io_utils import pdb_to_atoms
    masks, reports = [], []
    identity = None
    for path in paths:
        lines = path.read_text().splitlines()
        if sum(line.startswith('MODEL ') for line in lines) > 1:
            raise ValueError(f'{path}: use one PDB model per conformation.')
        rows = [line for line in lines if line.startswith(('ATOM  ', 'HETATM'))]
        if any(line[16:17].strip() for line in rows):
            raise ValueError(f'{path}: resolve alternate locations before splitting H.')
        if any(line[17:20].strip() in {'HOH', 'WAT', 'SOL', 'TIP', 'DOD'} for line in rows):
            raise ValueError(f'{path}: remove explicit water before solute SASA calculation.')
        atoms = pdb_to_atoms(path)
        if len(rows) != len(atoms):
            raise ValueError(f'{path}: malformed coordinate records.')
        current = [(line[12:27], atom[0]) for line, atom in zip(rows, atoms)]
        if identity is not None and current != identity:
            raise ValueError(f'{path}: atom identities/order differ between conformations.')
        identity = current
        mask, records = classify_hydrogens(atoms, **options)
        for record in records:
            line = rows[record['atom_index_1based'] - 1]
            record.update(pdb_serial=line[6:11].strip(), atom_name=line[12:16].strip(),
                          residue=line[17:20].strip(), chain=line[21:22].strip(),
                          residue_id=line[22:27].strip())
        masks.append(mask)
        reports.append({'pdb': str(path), 'hydrogens': records})
    return masks, reports
