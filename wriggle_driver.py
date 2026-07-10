#!/usr/bin/env python3
"""
wriggle_driver.py — gentle backbone bend + optional side-chain cleanup

This script launches wriggling jobs using stochastic_wriggle.py.

It keeps a single executable workflow, but stores recommended parameter
presets for different classes of systems.

Recommended presets
-------------------

1) Alpha-helical peptide
   Good starting values for short, mostly alpha-helical peptides:
       BB_SD = 6
       BB_WINDOW = 2
       BB_ATTEMPTS = 100
       BB_STRIDE = 2
       BB_CLASH_SCALE = 0.72
       CLASH_MODE = "on"
       SC_SD_BACKBONE = 5
       RUN_SC_OPT = True
       SC_OPT_PASSES = 6
       SC_OPT_SCALE = 0.98
       SC_OPT_VAR = 2
       USE_NO_SC_FLAG = False

2) Globular protein
   Good starting values for folded globular proteins:
       BB_SD = 5
       BB_WINDOW = 3
       BB_ATTEMPTS = 800
       BB_STRIDE = 4
       BB_CLASH_SCALE = 0.70
       CLASH_MODE = "off"
       SC_SD_BACKBONE = 0
       RUN_SC_OPT = True
       SC_OPT_PASSES = 4
       SC_OPT_SCALE = 0.98
       SC_OPT_VAR = 0.2
       USE_NO_SC_FLAG = True

Use:
    python wriggle_driver.py --input myprotein.pdb -n 20 --mode protein
    python wriggle_driver.py --input mypeptide.pdb -n 20 --mode peptide
"""

import argparse, os, sys, subprocess, pathlib
import json
# -------------------------------------------------------------------
# Presets
# -------------------------------------------------------------------

# Good starting values for mostly alpha-helical peptides
PEPTIDE_PRESET = {
    "BB_SD": 6,
    "BB_WINDOW": 2,
    "BB_ATTEMPTS": 100,
    "BB_STRIDE": 2,
    "BB_CLASH_SCALE": 0.72,
    "AXIS_REALIGN_FLAG": "--no-axis-realign",
    "CLASH_MODE": "on",
    "SCRIPT": "stochastic_wriggle.py",
    "SC_SD_BACKBONE": 5,
    "RUN_SC_OPT": True,
    "SC_OPT_PASSES": 6,
    "SC_OPT_SCALE": 0.98,
    "SC_OPT_VAR": 2,
    "USE_NO_SC_FLAG": False,
}

# Good starting values for globular proteins like colicin
PROTEIN_PRESET = {
    "BB_SD": 5,
    "BB_WINDOW": 3,
    "BB_ATTEMPTS": 800,
    "BB_STRIDE": 4,
    "BB_CLASH_SCALE": 0.70,
    "AXIS_REALIGN_FLAG": "--no-axis-realign",
    "CLASH_MODE": "off",
    "SCRIPT": "stochastic_wriggle.py",
    "SC_SD_BACKBONE": 0,
    "RUN_SC_OPT": True,
    "SC_OPT_PASSES": 4,
    "SC_OPT_SCALE": 0.98,
    "SC_OPT_VAR": 0.2,
    "USE_NO_SC_FLAG": True,
}


# no more hardcoded parameters, just some defaults
def main():
    arg_parse = argparse.ArgumentParser()
    arg_parse.add_argument("--input", default="peptide.pdb", help="Input PDB filename")
    arg_parse.add_argument("-n", type=int, default=20, help="Number of wriggled structures to generate")
    arg_parse.add_argument("--mode", choices=["peptide", "protein"], default="protein",
                           help="Select recommended preset parameters")

    # Optional overrides: start from preset, then change only selected parameters
    arg_parse.add_argument("--bb-sd", type=float, default=None)
    arg_parse.add_argument("--bb-window", type=int, default=None)
    arg_parse.add_argument("--bb-attempts", type=int, default=None)
    arg_parse.add_argument("--bb-stride", type=int, default=None)
    arg_parse.add_argument("--bb-clash-scale", type=float, default=None)

    arg_parse.add_argument("--sc-sd-backbone", type=float, default=None)
    arg_parse.add_argument("--sc-opt-passes", type=int, default=None)
    arg_parse.add_argument("--sc-opt-scale", type=float, default=None)
    arg_parse.add_argument("--sc-opt-var", type=float, default=None)

    args = arg_parse.parse_args()

    pdb_in = pathlib.Path(args.input)
    print("input=", pdb_in.name, "n=", args.n, "mode=", args.mode)

    # --- Choose preset without it getting changed by CLI---
    if args.mode == "peptide":
        preset = PEPTIDE_PRESET.copy()
    else:
        preset = PROTEIN_PRESET.copy()

# --- Optional command-line overrides ---
# These allow users to start from a preset and modify only selected parameters.

    if args.bb_sd is not None:
        preset["BB_SD"] = args.bb_sd

    if args.bb_window is not None:
        preset["BB_WINDOW"] = args.bb_window

    if args.bb_attempts is not None:
        preset["BB_ATTEMPTS"] = args.bb_attempts

    if args.bb_stride is not None:
        preset["BB_STRIDE"] = args.bb_stride

    if args.bb_clash_scale is not None:
        preset["BB_CLASH_SCALE"] = args.bb_clash_scale

    if args.sc_sd_backbone is not None:
        preset["SC_SD_BACKBONE"] = args.sc_sd_backbone

    if args.sc_opt_passes is not None:
        preset["SC_OPT_PASSES"] = args.sc_opt_passes

    if args.sc_opt_scale is not None:
        preset["SC_OPT_SCALE"] = args.sc_opt_scale

    if args.sc_opt_var is not None:
        preset["SC_OPT_VAR"] = args.sc_opt_var

    with open("wriggle_parameters_used.json", "w") as f:
        json.dump(
            {
                "mode": args.mode,
                "n_structures": args.n,
                "parameters": preset,
            },
            f,
            indent=2,
        )


    # Backbone
    BB_SD = preset["BB_SD"]
    BB_WINDOW = preset["BB_WINDOW"]
    BB_ATTEMPTS = preset["BB_ATTEMPTS"]
    BB_STRIDE = preset["BB_STRIDE"]
    BB_CLASH_SCALE = preset["BB_CLASH_SCALE"]
    AXIS_REALIGN_FLAG = preset["AXIS_REALIGN_FLAG"]

    # Clash mode
    CLASH_MODE = preset["CLASH_MODE"]
    SCRIPT = preset["SCRIPT"]

    # --- Move side chains slightly during backbone pass ---
    SC_SD_BACKBONE = preset["SC_SD_BACKBONE"]

    # --- Side chain optimiser ---
    RUN_SC_OPT = preset["RUN_SC_OPT"]
    SC_OPT_PASSES = preset["SC_OPT_PASSES"]
    SC_OPT_SCALE = preset["SC_OPT_SCALE"]
    SC_OPT_VAR = preset["SC_OPT_VAR"]

    # --- Whether to skip initial side-chain sampling during backbone step ---
    USE_NO_SC_FLAG = preset["USE_NO_SC_FLAG"]

    # --- Omega tweaks: off in this workflow ---
    OMEGA_FLAGS = ["--no-omega"]

    print(f"[launch] backbone: sd={BB_SD} window={BB_WINDOW} attempts={BB_ATTEMPTS} "
          f"stride={BB_STRIDE} clash_scale={BB_CLASH_SCALE}")
    print(f"[launch] clash mode: {CLASH_MODE}")
    print(f"[launch] initial side-chain sampling skipped: {USE_NO_SC_FLAG}")
    if RUN_SC_OPT:
        print(f"[launch] sc-opt: passes={SC_OPT_PASSES} scale={SC_OPT_SCALE} variation={SC_OPT_VAR}")
    else:
        print("[launch] sc-opt: disabled")

    tmp_dir = pathlib.Path("tmp_launch")
    tmp_dir.mkdir(exist_ok=True)

    out_dir = pathlib.Path("structures")
    out_dir.mkdir(exist_ok=True)



    for i in range(1, args.n + 1):
        print(f"[launch][{i}/{args.n}] backbone bend")
        tmp_bb = tmp_dir / f"tmp_bb_{i}.pdb"
        out_pdb = out_dir / f"conf_{i:03d}.pdb"

        cmd_bb = [
            sys.executable, SCRIPT, pdb_in.name,
            "-o", str(tmp_bb),
            "--no-omega",
            "--sc-sd", str(SC_SD_BACKBONE),
            "--bb-sd", str(BB_SD),
            "--bb-window", str(BB_WINDOW),
            "--bb-attempts", str(BB_ATTEMPTS),
            "--bb-stride", str(BB_STRIDE),
            "--clash-mode", str(CLASH_MODE),
            "--bb-clash-scale", str(BB_CLASH_SCALE),
            AXIS_REALIGN_FLAG,
            "--report",
            "--seed", str(i),
        ]

        # For globular proteins we usually skip the initial side-chain rotamer shake
        if USE_NO_SC_FLAG:
            cmd_bb.append("--no-sc")

        subprocess.run(cmd_bb, check=True)

        cmd_sc = [
            sys.executable, SCRIPT, str(tmp_bb),
            "-o", str(out_pdb),
            "--no-backbone", "--no-omega",
            "--sc-sd", "0",
            "--sc-optimise",
            "--sc-opt-passes", str(SC_OPT_PASSES),
            "--sc-opt-scale", str(SC_OPT_SCALE),
            "--sc-opt-variation", str(SC_OPT_VAR),
            "--report",
            "--seed", str(i + 1000)
        ]

        # If the backbone stage used --no-sc, keep the next stage consistent
        if USE_NO_SC_FLAG:
            cmd_sc.append("--no-sc")

        subprocess.run(cmd_sc, check=True)

if __name__ == "__main__":
    main()
