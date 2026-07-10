#!/usr/bin/env python3
"""
io_utils.py

Helper functions for reading and writing simple coordinate formats used in the
NR workflow.

This file collects only the I/O functions actually needed by the current
pipeline.
"""

from pathlib import Path


# ------------------------ XYZ helpers ------------------------

def write_xyz(path, atoms):
    """
    Write an XYZ file.

    Parameters
    ----------
    path : str or Path
        Output XYZ filename.

    atoms : list of tuples
        [(elem, x, y, z), ...]
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    with open(path, "w") as f:
        f.write(str(len(atoms)) + "\n")
        f.write(str(path) + "\n")
        for elem, x, y, z in atoms:
            f.write(f"{elem} {x:.5f} {y:.5f} {z:.5f}\n")


def write_xyz_trajectory(path, frames):
    """
    Write a multi-frame XYZ trajectory.

    frames : list
        [
            [(elem, x, y, z), ...],   # frame 1
            [(elem, x, y, z), ...],   # frame 2
            ...
        ]

    VMD should read this as a trajectory with one frame per conformation.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    with open(path, "w") as f:
        for frame_index, atoms in enumerate(frames):
            f.write(str(len(atoms)) + "\n")
            f.write(f"frame {frame_index}\n")

            for elem, x, y, z in atoms:
                f.write(f"{elem} {x:.5f} {y:.5f} {z:.5f}\n")

# ------------------------ PDB / element helpers ------------------------

def infer_element_from_pdb_line(line):
    """
    Infer element symbol from a PDB ATOM/HETATM line.

    First tries the standard element column (77–78).
    If empty, falls back to the first alphabetic character of the atom name.
    """
    if len(line) >= 78:
        elem = line[76:78].strip()
        if elem:
            return elem[0].upper()

    name = line[12:16].strip()
    for ch in name:
        if ch.isalpha():
            return ch.upper()

    return "C"


def pdb_to_xyz(pdb_path, xyz_path):
    """
    Convert a PDB file to XYZ format.

    Parameters
    ----------
    pdb_path : str or Path
        Input PDB filename.

    xyz_path : str or Path
        Output XYZ filename.
    """
    pdb_path = Path(pdb_path)
    xyz_path = Path(xyz_path)

    with open(pdb_path) as f:
        lines = f.readlines()

    out = []
    for line in lines:
        if line.startswith(("ATOM  ", "HETATM")) and len(line) >= 54:
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except ValueError:
                continue

            elem = infer_element_from_pdb_line(line)
            out.append(f"{elem} {x:.5f} {y:.5f} {z:.5f}\n")

    xyz_path.parent.mkdir(parents=True, exist_ok=True)
    with open(xyz_path, "w") as f:
        f.writelines(out)


def pdb_to_atoms(pdb_path):
    """
    Read atomic coordinates from a PDB file.

    Returns
    -------
    atoms : list of tuples
        [(elem, x, y, z), ...]
    """
    pdb_path = Path(pdb_path)

    atoms = []
    with open(pdb_path) as f:
        for line in f:
            if line.startswith(("ATOM  ", "HETATM")) and len(line) >= 54:
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                except ValueError:
                    continue

                elem = infer_element_from_pdb_line(line)
                atoms.append((elem, x, y, z))

    return atoms


def write_pdb_from_template(template_pdb, out_pdb, atoms):
    """
    Write a PDB file using a template PDB for all fields except coordinates.

    Parameters
    ----------
    template_pdb : str or Path
        Path to the original PDB file. This file provides:
        - atom names
        - residue names and IDs
        - chain IDs
        - occupancies, temp factors, element symbols, etc.

    out_pdb : str or Path
        Path where the new PDB will be written.

    atoms : list of tuples
        [(elem, x, y, z), ...]

        IMPORTANT:
        - Must be in the same atom order as template_pdb
        - Only coordinates are replaced
    """
    template_pdb = Path(template_pdb)
    out_pdb = Path(out_pdb)

    with open(template_pdb, "r") as fin, open(out_pdb, "w") as fout:
        atom_index = 0

        for line in fin:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                atom = atoms[atom_index]

                if len(atom) == 4:
                    _, x, y, z = atom
                elif len(atom) == 3:
                    x, y, z = atom
                else:
                    raise ValueError("Atom format must be (x,y,z) or (elem,x,y,z)")

                newline = (
                    line[:30] +
                    f"{x:8.3f}{y:8.3f}{z:8.3f}" +
                    line[54:]
                )

                fout.write(newline)
                atom_index += 1
            else:
                fout.write(line)

    if atom_index != len(atoms):
        raise ValueError(
            f"Number of atoms written ({atom_index}) "
            f"does not match provided atoms ({len(atoms)})"
        )
