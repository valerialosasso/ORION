"""PQR -> standard fixed-column PDB conversion and ORION input validation."""
from pathlib import Path
import math


def _element(atom_name):
    # PDB2PQR protein atom names: H*, C*, N*, O*, S*, P*.
    # For unusual ions/ligands, require a dedicated chemical parser.
    name = atom_name.lstrip('0123456789').upper()
    if not name or name[0] not in 'HCNOSP':
        raise ValueError(f'Cannot safely infer element for PQR atom {atom_name!r}.')
    return name[0]


def convert_pqr_to_pdb(source, destination):
    """Convert common PDB2PQR whitespace records (with/without chain) to PDB.

    Charges/radii are deliberately not transferred; ORION does not use them.
    """
    source, destination = Path(source), Path(destination)
    rows = []
    for line_no, line in enumerate(source.read_text().splitlines(), 1):
        if not line.startswith(('ATOM  ', 'HETATM')):
            continue
        fields = line.split()
        if len(fields) == 10:
            rec, serial, atom, res, resid, x, y, z, charge, radius = fields
            chain = 'A'
        elif len(fields) == 11:
            rec, serial, atom, res, chain, resid, x, y, z, charge, radius = fields
        else:
            raise ValueError(f'{source}:{line_no}: unsupported PQR record ({len(fields)} fields).')
        try:
            serial, resid = int(serial), int(resid)
            x, y, z = float(x), float(y), float(z)
            float(charge), float(radius)
        except ValueError as exc:
            raise ValueError(f'{source}:{line_no}: invalid PQR numerical field') from exc
        if not all(map(math.isfinite, (x,y,z))):
            raise ValueError(f'{source}:{line_no}: nonfinite coordinates')
        if not (1 <= serial <= 99999 and -999 <= resid <= 9999):
            raise ValueError(f'{source}:{line_no}: PDB serial/residue number out of range')
        if len(atom)>4 or len(res)>3 or len(chain)!=1:
            raise ValueError(f'{source}:{line_no}: PDB identifier cannot be represented')
        element = _element(atom)
        atom_field = f' {atom:<3}' if len(atom)<4 else atom
        rows.append(f'{rec:<6}{serial:5d} {atom_field:4s} {res:>3s} {chain}{resid:4d}    '
                    f'{x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{0.0:6.2f}'
                    f'          {element:>2s}\n')
    if not rows:
        raise ValueError(f'{source}: no PQR ATOM/HETATM records found')
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text('REMARK Converted from PQR; charge and radius omitted\n' + ''.join(rows) + 'END\n')
    return destination


def validate_pdb(path, require_hydrogens=False):
    from io_utils import pdb_to_atoms
    path = Path(path)
    lines = [line for line in path.read_text().splitlines() if line.startswith(('ATOM  ', 'HETATM'))]
    if not lines:
        raise ValueError(f'{path}: no coordinate records')
    if any(line[16:17].strip() for line in lines):
        raise ValueError(f'{path}: resolve alternate locations before wriggling.')
    atoms = pdb_to_atoms(path)
    if len(atoms) != len(lines):
        raise ValueError(f'{path}: malformed PDB coordinate records')
    if require_hydrogens and not any(atom[0] == 'H' for atom in atoms):
        raise ValueError(f'{path}: no explicit hydrogens. Add hydrogens with PDB2PQR or another preparation tool before using --split-hydrogens.')
    return len(atoms)
