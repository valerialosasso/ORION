#!/usr/bin/env bash
set -euo pipefail

# Minimal end-to-end ORION test.
# Generates one conformation and one orientation, then checks key outputs.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

PYTHON_BIN="${PYTHON_BIN:-python3}"
INPUT_PDB="${REPO_ROOT}/examples/Colicin-N_1A87.pdb"
OUTPUT_DIR="${REPO_ROOT}/test_smoke"
LOOKUP_FILE="${OUTPUT_DIR}/smoke_lookup.json"

if [[ ! -f "${REPO_ROOT}/run_pipeline.py" ]]; then
    echo "ERROR: run_pipeline.py was not found in ${REPO_ROOT}" >&2
    exit 1
fi

if [[ ! -f "${INPUT_PDB}" ]]; then
    echo "ERROR: smoke-test PDB was not found at ${INPUT_PDB}" >&2
    exit 1
fi

rm -rf "${OUTPUT_DIR}"

echo "Running ORION smoke test..."

"${PYTHON_BIN}" "${REPO_ROOT}/run_pipeline.py" \
    --pdb "${INPUT_PDB}" \
    --theta-angles 0:0:1 \
    --phi-angles 0:0:1 \
    --outdir "${OUTPUT_DIR}" \
    --json smoke_lookup.json \
    --box-size 150 \
    --bin 1.5 \
    --nstruc 1 \
    --wriggle-mode protein

required_files=(
    "${LOOKUP_FILE}"
    "${OUTPUT_DIR}/wriggles_base/structures/conf_001.pdb"
    "${OUTPUT_DIR}/theta0_phi0/tempFiles/trajectory.xyz"
    "${OUTPUT_DIR}/theta0_phi0/density_C.dat"
    "${OUTPUT_DIR}/theta0_phi0/density_H.dat"
    "${OUTPUT_DIR}/theta0_phi0/density_N.dat"
    "${OUTPUT_DIR}/theta0_phi0/density_O.dat"
    "${OUTPUT_DIR}/theta0_phi0/vol_frac_along_z.dat"
)

for file in "${required_files[@]}"; do
    if [[ ! -s "${file}" ]]; then
        echo "SMOKE TEST FAILED: missing or empty file: ${file}" >&2
        exit 1
    fi
done

"${PYTHON_BIN}" - "${LOOKUP_FILE}" <<'PY'
import json
import sys
from pathlib import Path

path = Path(sys.argv[1])
with path.open() as handle:
    data = json.load(handle)

if not isinstance(data, list) or not data:
    raise SystemExit("SMOKE TEST FAILED: lookup JSON is empty or has an unexpected format")

print("Lookup JSON is valid.")
PY

echo "SMOKE TEST PASSED"
