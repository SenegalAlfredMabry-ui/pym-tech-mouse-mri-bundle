#!/bin/bash -l
#SBATCH -J unzip_mouse_mri
#SBATCH -p cloud
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --time=08:00:00
#SBATCH -o sbatch-unzip-mouse-mri-%j.out
#SBATCH -e sbatch-unzip-mouse-mri-%j.err

set -euo pipefail

REPO_ROOT=${SNHEART_ROOT:-/path/to/SNHEART}
ZIP_PATH="${REPO_ROOT}/MRI_data/Raw/nonhuman_mouse/Mouse_MRI_Directory.zip"
OUT_ROOT="${REPO_ROOT}/MRI_data/Raw/nonhuman_mouse"
EXTRACT_DIR="${OUT_ROOT}/Mouse_MRI_Directory"
REPORT_DIR="${OUT_ROOT}/unzip_reports"
MARKER="${OUT_ROOT}/Mouse_MRI_Directory.UNZIP_COMPLETE.txt"

mkdir -p "${REPORT_DIR}"
REPORT_PATH="${REPORT_DIR}/unzip_mouse_mri_directory_${SLURM_JOB_ID}.txt"

echo "[$(date -u '+%F %T UTC')] Starting unzip job ${SLURM_JOB_ID}" | tee "${REPORT_PATH}"
echo "ZIP_PATH=${ZIP_PATH}" | tee -a "${REPORT_PATH}"
echo "OUT_ROOT=${OUT_ROOT}" | tee -a "${REPORT_PATH}"

if [[ ! -f "${ZIP_PATH}" ]]; then
  echo "ERROR: zip file not found: ${ZIP_PATH}" | tee -a "${REPORT_PATH}" >&2
  exit 1
fi

cd "${OUT_ROOT}"
export UNZIP_DISABLE_ZIPBOMB_DETECTION=TRUE
unzip -q -o "${ZIP_PATH}" -d "${OUT_ROOT}"

ZIP_LIST="${REPORT_DIR}/zip_manifest_${SLURM_JOB_ID}.txt"
EXT_LIST="${REPORT_DIR}/extracted_manifest_${SLURM_JOB_ID}.txt"
MISSING_LIST="${REPORT_DIR}/missing_from_extract_${SLURM_JOB_ID}.txt"

unzip -Z -1 "${ZIP_PATH}" | sed '/\/$/d' | sort > "${ZIP_LIST}"
find "Mouse_MRI_Directory" -type f | sort > "${EXT_LIST}"
comm -23 "${ZIP_LIST}" "${EXT_LIST}" > "${MISSING_LIST}"

zip_count=$(wc -l < "${ZIP_LIST}" | tr -d ' ')
ext_count=$(wc -l < "${EXT_LIST}" | tr -d ' ')
missing_count=$(wc -l < "${MISSING_LIST}" | tr -d ' ')

echo "zip_file_count=${zip_count}" | tee -a "${REPORT_PATH}"
echo "extracted_file_count=${ext_count}" | tee -a "${REPORT_PATH}"
echo "missing_file_count=${missing_count}" | tee -a "${REPORT_PATH}"

if [[ "${missing_count}" != "0" ]]; then
  echo "ERROR: extraction incomplete; missing files listed at ${MISSING_LIST}" | tee -a "${REPORT_PATH}" >&2
  exit 2
fi

{
  echo "completed_utc=$(date -u '+%F %T UTC')"
  echo "job_id=${SLURM_JOB_ID}"
  echo "zip_file_count=${zip_count}"
  echo "extracted_file_count=${ext_count}"
  echo "report_path=${REPORT_PATH}"
} > "${MARKER}"

echo "SUCCESS: extraction complete and verified" | tee -a "${REPORT_PATH}"
echo "Marker written: ${MARKER}" | tee -a "${REPORT_PATH}"
