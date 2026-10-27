"""PDB2PQR -> PDB conversion for ORION, including four-character residue names.

The PQR input remains unchanged. Explicit mappings preserve coordinates and atom
order; patch/protonation metadata is reported but cannot be represented in PDB.
"""
from pathlib import Path
import math
from collections import Counter

RESIDUE_MAP = {'DISU': 'CYS', 'GLUP': 'GLU'}


def _element(atom_name):
    name = atom_name.lstrip('0123456789').upper()
    if not name or name[0] not in 'HCNOSP':
        raise ValueError(f'Cannot safely infer element for PQR atom {atom_name!r}.')
    return name[0]


def convert_pqr_to_pdb(source, destination):
    source, destination = Path(source), Path(destination)
    rows, mappings = [], Counter()
    for line_no, line in enumerate(source.read_text().splitlines(), 1):
        if not line.startswith(('ATOM  ', 'HETATM')):
            continue
        fields = line.split()
        if len(fields) in (10, 11):
            if len(fields) == 10:
                rec, serial, atom, res, resid, x, y, z, charge, radius = fields
                chain = 'A'
            else:
                rec, serial, atom, res, chain, resid, x, y, z, charge, radius = fields
        elif len(fields) == 9:
            # Some PDB2PQR variants write a four-character residue name in
            # columns 17-20, with no separating space: 1CBDISU / OE2GLUP.
            # Parse the fixed-position identifier fields, NOT the fused token.
            rec = line[:6].strip()
            serial = line[6:11].strip()
            atom = line[12:16].strip()
            res = line[16:20].strip() if line[16] != ' ' else line[17:21].strip()
            chain = line[21].strip() or 'A'
            resid = line[22:26].strip()
            x, y, z, charge, radius = fields[-5:]
            if res not in RESIDUE_MAP:
                raise ValueError(f'{source}:{line_no}: unrecognised four-character residue {res!r}; review manually')
        else:
            raise ValueError(f'{source}:{line_no}: unsupported PQR record ({len(fields)} fields).')
        try:
            serial, resid = int(serial), int(resid)
            x, y, z = float(x), float(y), float(z)
            charge, radius = float(charge), float(radius)
        except ValueError as exc:
            raise ValueError(f'{source}:{line_no}: invalid numerical field') from exc
        if not all(map(math.isfinite, (x, y, z, charge, radius))):
            raise ValueError(f'{source}:{line_no}: nonfinite PQR numerical field')
        if not (1 <= serial <= 99999 and -999 <= resid <= 9999):
            raise ValueError(f'{source}:{line_no}: PDB serial/residue number out of range')
        if res in RESIDUE_MAP:
            mappings[(res, RESIDUE_MAP[res])] += 1
            res = RESIDUE_MAP[res]
        if atom in ('1CB', '1SG') and res == 'CYS':
            mappings[(atom, atom[1:])] += 1
            atom = atom[1:]
        if len(atom) > 4 or len(res) > 3 or len(chain) != 1:
            raise ValueError(f'{source}:{line_no}: identifier cannot be represented in PDB')
        element = _element(atom)
        atom_field = f' {atom:<3}' if len(atom) < 4 else atom
        rows.append(f'{rec:<6}{serial:5d} {atom_field:4s} {res:>3s} {chain}{resid:4d}    '
                    f'{x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{0.0:6.2f}'
                    f'          {element:>2s}\n')
    if not rows:
        raise ValueError(f'{source}: no PQR coordinate records')
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text('REMARK Converted from PQR; charges/radii omitted\n' + ''.join(rows) + 'END\n')
    if mappings:
        print('PQR residue/atom-name normalization:', dict(mappings))
    return destination


def validate_pdb(path, require_hydrogens=False):
    from io_utils import pdb_to_atoms
    path = Path(path)
    lines = [line for line in path.read_text().splitlines() if line.startswith(('ATOM  ', 'HETATM'))]
    if not lines:
        raise ValueError(f'{path}: no coordinate records')
    if any(line[16:17].strip() for line in lines):
        raise ValueError(f'{path}: resolve alternate locations before wriggling')
    atoms = pdb_to_atoms(path)
    if len(atoms) != len(lines):
        raise ValueError(f'{path}: malformed PDB coordinate records')
    if require_hydrogens and not any(atom[0] == 'H' for atom in atoms):
        raise ValueError(f'{path}: no explicit hydrogens. Add them before --split-hydrogens.')
    return len(atoms)
