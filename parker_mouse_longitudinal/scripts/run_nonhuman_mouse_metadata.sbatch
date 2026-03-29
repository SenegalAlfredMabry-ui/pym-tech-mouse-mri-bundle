#!/bin/bash -l
#SBATCH -J nonhuman_mouse_meta
#SBATCH -p cloud
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --time=01:00:00
#SBATCH --output=${SNHEART_ROOT:-/path/to/SNHEART}/Logs/sbatch-nonhuman-metadata-%j.out
#SBATCH --error=${SNHEART_ROOT:-/path/to/SNHEART}/Logs/sbatch-nonhuman-metadata-%j.err

set -euo pipefail

repo=${SNHEART_ROOT:-/path/to/SNHEART}
root=${NONHUMAN_MOUSE_ROOT:-$repo/MRI_data/Raw/nonhuman_mouse}

echo "[INFO] Start nonhuman mouse metadata extraction: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "[INFO] root=$root"

"$repo/venv/bin/python" "$repo/Scripts/extract_nonhuman_mouse_dicom_metadata.py" \
  --root "$root"

echo "[INFO] End nonhuman mouse metadata extraction: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
