#!/bin/bash -l
#SBATCH -J mouse_dcm2nii
#SBATCH -p cloud
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --time=24:00:00
#SBATCH -o sbatch-mouse-dcm2nii-%j.out
#SBATCH -e sbatch-mouse-dcm2nii-%j.err

set -euo pipefail

REPO_ROOT=${SNHEART_ROOT:-/path/to/SNHEART}
cd "${REPO_ROOT}"

MANIFEST_TSV=${MANIFEST_TSV:-${REPO_ROOT}/MRI_data/Raw/nonhuman_mouse/Mouse_MRI_Directory_organized/manifests/organized_manifest.tsv}
OUT_ROOT=${OUT_ROOT:-${REPO_ROOT}/MRI_data/Raw/nonhuman_mouse/NIFTI_longitudinal}
DCM2NIIX_PATH=${DCM2NIIX_PATH:-${DCM2NIIX_PATH:-dcm2niix}}
INCLUDE_MODALITIES=${INCLUDE_MODALITIES:-T1gre,T2_RARE,T2tse,fMRI_6000,fMRI_7000,TOF}
OVERWRITE=${OVERWRITE:-0}
DRY_RUN=${DRY_RUN:-0}

cmd=(
  python Scripts/convert_nonhuman_mouse_organized_to_nifti.py
  --manifest-tsv "${MANIFEST_TSV}"
  --out-root "${OUT_ROOT}"
  --dcm2niix-path "${DCM2NIIX_PATH}"
  --include-modalities "${INCLUDE_MODALITIES}"
)

if [[ "${OVERWRITE}" == "1" ]]; then
  cmd+=(--overwrite)
fi
if [[ "${DRY_RUN}" == "1" ]]; then
  cmd+=(--dry-run)
fi

echo "Running command:"
printf '%q ' "${cmd[@]}"
echo
"${cmd[@]}"
