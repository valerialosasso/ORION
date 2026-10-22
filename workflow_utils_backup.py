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
    "density_P.dat", "density_H_group1.dat", "density_H_group2.dat",
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
    """Make the wriggling script available, repairing stale or broken links."""
    orient_dir = Path(orient_dir)
    orient_dir.mkdir(parents=True, exist_ok=True)
    stoch_src = (Path(__file__).resolve().parent / script_name).resolve(strict=True)
    stoch_dst = orient_dir / script_name

    if stoch_dst.is_symlink():
        if stoch_dst.exists() and stoch_dst.resolve() == stoch_src:
            return
        stoch_dst.unlink()  # remove stale/broken link, not its target
    elif stoch_dst.exists():
        if stoch_dst.is_dir():
            raise IsADirectoryError(f"Expected a file at {stoch_dst}")
        # Preserve a pre-existing regular file, as the original function did.
        return

    try:
        stoch_dst.symlink_to(stoch_src)
    except OSError:
        # Only fall back when no destination was created.
        if stoch_dst.is_symlink() or stoch_dst.exists():
            raise
        shutil.copyfile(stoch_src, stoch_dst)
