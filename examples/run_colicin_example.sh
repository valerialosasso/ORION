#!/usr/bin/env bash
set -euo pipefail

# Run the full Colicin-N example from anywhere inside the repository.
# The script resolves paths relative to its own location.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

PYTHON_BIN="${PYTHON_BIN:-python3}"
INPUT_PDB="${REPO_ROOT}/examples/Colicin-N_1A87.pdb"
OUTPUT_DIR="${REPO_ROOT}/example_output/colicin_orientations"

if [[ ! -f "${REPO_ROOT}/run_pipeline.py" ]]; then
    echo "ERROR: run_pipeline.py was not found in ${REPO_ROOT}" >&2
    exit 1
fi

if [[ ! -f "${INPUT_PDB}" ]]; then
    echo "ERROR: example PDB was not found at ${INPUT_PDB}" >&2
    exit 1
fi

echo "Running ORION Colicin-N example"
echo "Python: ${PYTHON_BIN}"
echo "Input:  ${INPUT_PDB}"
echo "Output: ${OUTPUT_DIR}"

"${PYTHON_BIN}" "${REPO_ROOT}/run_pipeline.py" \
    --pdb "${INPUT_PDB}" \
    --theta-angles 0:90:5 \
    --phi-angles 0:90:5 \
    --outdir "${OUTPUT_DIR}" \
    --json colicin_N_lookup.json \
    --box-size 200 \
    --bin 1.5 \
    --nstruc 100 \
    --wriggle-mode protein

echo "ORION example completed successfully."
echo "Lookup table: ${OUTPUT_DIR}/colicin_N_lookup.json"
