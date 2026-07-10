#!/usr/bin/env python3
"""
workflow_utils.py

Helper functions for cleaning, angle parsing, and workflow setup.
"""

import shutil
from pathlib import Path


# ------------------------ constants ---------------------------------

# files that get produced during density / SLD calculation
DISTRIBUTION_FILES = [
    "density_C.dat", "density_H.dat", "density_N.dat", "density_O.dat",
    "density_P.dat",
    "vol_frac_along_z.dat",
    "sld_C.dat", "sld_H.dat", "sld_N.dat", "sld_O.dat", "sld_total.dat"
]


# ------------------------ cleaning ---------------------------------

def clean_distributions(root):
    """
    Remove old output distributions and temporary files.

    This is used to avoid mixing outputs from previous runs with the current one.
    """
    root = Path(root)

    # remove density / SLD files if they exist
    for file_name in DISTRIBUTION_FILES:
        file_path = root / file_name
        if file_path.exists():
            try:
                file_path.unlink()
            except Exception:
                pass

    # remove temporary directory if present
    temp_dir = root / "tempFiles"
    if temp_dir.exists():
        try:
            shutil.rmtree(temp_dir)
        except Exception:
            pass


# ------------------------ angle parsing ------------------------

def parse_angle_spec(spec):
    """
    Format = 'a:b:c' (inclusive range with step)

    Example:
        "0:90:5" → [0, 5, 10, ..., 90]

    Returns list of integers (angles)
    """
    spec = str(spec).strip()

    a, b, c = spec.split(":")
    a = float(a)
    b = float(b)
    c = float(c)

    values = []
    v = a

    # include endpoint if exactly reached
    while v <= b:
        values.append(int(round(v)))
        v += c

    return values


# ------------------------ stochastic script handling ------------------------

def ensure_stochastic_in_dir(orient_dir, script_name="stochastic_wriggle.py"):
    """
    Ensure stochastic_wriggle.py is available inside orient_dir.

    This is needed because the wriggling step is executed inside each
    orientation directory.

    Strategy:
      1. Try to create a symlink (cleaner, avoids duplication)
      2. If that fails (e.g. filesystem restrictions), copy the file
    """
    orient_dir = Path(orient_dir)

    # locate script in current code directory
    script_dir = Path(__file__).resolve().parent
    stoch_src = script_dir / script_name
    stoch_dst = orient_dir / script_name

    # only create if missing
    if not stoch_dst.exists():
        try:
            # try symlink first
            stoch_dst.symlink_to(stoch_src.resolve())
        except Exception:
            # fallback: copy file
            shutil.copyfile(stoch_src.resolve(), stoch_dst)
