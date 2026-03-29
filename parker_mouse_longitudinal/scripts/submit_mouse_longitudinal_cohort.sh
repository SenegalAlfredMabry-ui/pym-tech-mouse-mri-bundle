#!/bin/bash -l
# Submit one longitudinal mouse job per cohort row.
#
# Cohort TSV format (tab-separated):
#   mouse_id<TAB>mouse_root<TAB>timepoint_glob(optional)
# Header row "mouse_id<TAB>mouse_root..." is allowed.

set -euo pipefail

usage() {
  cat <<'EOF'
Usage: Scripts/submit_mouse_longitudinal_cohort.sh --cohort-tsv <path> [options]

Required:
  --cohort-tsv <path>       TSV with: mouse_id, mouse_root, optional timepoint_glob

Optional run settings:
  --methods <csv>           Default: ants,sswarper
  --no-parcels <0|1>        Default: 1
  --atlas-parcels <path>    Required if --no-parcels 0
  --atlas-template <path>   Optional; pipeline default is used when omitted
  --atlas-labels <path>     Optional labels table
  --out-root <dir>          Default: MRI_data/processed/nonhuman_longitudinal
  --tp-glob <glob>          Default timepoint glob when row does not provide one (default: *)
  --threads <n>             Default: 4
  --template-n4 <0|1>       Default: 1
  --driver-prefix-priority <csv>   Default: 3000,9000,4000
  --driver-exclude-prefixes <csv>  Default: 6000,7000,8000

Optional sbatch settings:
  --partition <name>        Default: cloud
  --cpus <n>                Default: 4
  --time <HH:MM:SS>         Default: 48:00:00
  --dry-run                 Print submit commands without launching jobs

EOF
}

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(cd -- "${SCRIPT_DIR}/.." && pwd)

COHORT_TSV=""
METHODS="ants,sswarper"
NO_PARCELS=1
ATLAS_PARCELS=""
ATLAS_TEMPLATE=""
ATLAS_LABELS=""
OUT_ROOT="${REPO_ROOT}/MRI_data/processed/nonhuman_longitudinal"
TP_GLOB_DEFAULT="*"
THREADS=4
TEMPLATE_N4=1
DRIVER_PREFIX_PRIORITY="3000,9000,4000"
DRIVER_EXCLUDE_PREFIXES="6000,7000,8000"

PARTITION="cloud"
CPUS=4
TIME_LIMIT="48:00:00"
DRY_RUN=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --cohort-tsv) COHORT_TSV=$2; shift 2 ;;
    --methods) METHODS=$2; shift 2 ;;
    --no-parcels) NO_PARCELS=$2; shift 2 ;;
    --atlas-parcels) ATLAS_PARCELS=$2; shift 2 ;;
    --atlas-template) ATLAS_TEMPLATE=$2; shift 2 ;;
    --atlas-labels) ATLAS_LABELS=$2; shift 2 ;;
    --out-root) OUT_ROOT=$2; shift 2 ;;
    --tp-glob) TP_GLOB_DEFAULT=$2; shift 2 ;;
    --threads) THREADS=$2; shift 2 ;;
    --template-n4) TEMPLATE_N4=$2; shift 2 ;;
    --driver-prefix-priority) DRIVER_PREFIX_PRIORITY=$2; shift 2 ;;
    --driver-exclude-prefixes) DRIVER_EXCLUDE_PREFIXES=$2; shift 2 ;;
    --partition) PARTITION=$2; shift 2 ;;
    --cpus) CPUS=$2; shift 2 ;;
    --time) TIME_LIMIT=$2; shift 2 ;;
    --dry-run) DRY_RUN=1; shift ;;
    --help|-h) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage; exit 1 ;;
  esac
done

[[ -n ${COHORT_TSV} ]] || { echo "ERROR: --cohort-tsv is required" >&2; exit 1; }
COHORT_TSV=$(realpath "${COHORT_TSV}")
[[ -f ${COHORT_TSV} ]] || { echo "ERROR: cohort TSV not found: ${COHORT_TSV}" >&2; exit 1; }

if [[ ${NO_PARCELS} -eq 0 && -z ${ATLAS_PARCELS} ]]; then
  echo "ERROR: --atlas-parcels is required when --no-parcels 0" >&2
  exit 1
fi
if [[ -n ${ATLAS_PARCELS} ]]; then
  ATLAS_PARCELS=$(realpath "${ATLAS_PARCELS}")
fi
if [[ -n ${ATLAS_TEMPLATE} ]]; then
  ATLAS_TEMPLATE=$(realpath "${ATLAS_TEMPLATE}")
fi
if [[ -n ${ATLAS_LABELS} ]]; then
  ATLAS_LABELS=$(realpath "${ATLAS_LABELS}")
fi

STAMP=$(date +%Y%m%dT%H%M%S)
SUBMIT_DIR="${OUT_ROOT}/cohort_submissions"
mkdir -p "${SUBMIT_DIR}"
SUBMIT_TSV="${SUBMIT_DIR}/mouse_longitudinal_submit_${STAMP}.tsv"
printf 'mouse_id\tmouse_root\ttp_glob\tjob_id\tsubmit_state\n' > "${SUBMIT_TSV}"

echo "Cohort TSV: ${COHORT_TSV}"
echo "Submission log: ${SUBMIT_TSV}"

submit_count=0
skip_count=0

while IFS=$'\t' read -r mouse_id mouse_root tp_glob _rest; do
  [[ -z ${mouse_id} ]] && continue
  [[ ${mouse_id} =~ ^# ]] && continue
  if [[ ${mouse_id} == "mouse_id" ]]; then
    continue
  fi
  if [[ -z ${mouse_root} ]]; then
    echo "WARN: skipping row with empty mouse_root for ${mouse_id}" >&2
    printf '%s\t%s\t%s\t\tSKIPPED_EMPTY_ROOT\n' "${mouse_id}" "${mouse_root}" "${tp_glob:-${TP_GLOB_DEFAULT}}" >> "${SUBMIT_TSV}"
    skip_count=$((skip_count + 1))
    continue
  fi

  if [[ ! -d ${mouse_root} ]]; then
    echo "WARN: skipping ${mouse_id}; missing directory: ${mouse_root}" >&2
    printf '%s\t%s\t%s\t\tSKIPPED_MISSING_DIR\n' "${mouse_id}" "${mouse_root}" "${tp_glob:-${TP_GLOB_DEFAULT}}" >> "${SUBMIT_TSV}"
    skip_count=$((skip_count + 1))
    continue
  fi

  row_tp_glob=${tp_glob:-${TP_GLOB_DEFAULT}}
  env_items=(
    "MOUSE_ROOT=${mouse_root}"
    "MOUSE_ID=${mouse_id}"
    "METHODS=${METHODS}"
    "NO_PARCELS=${NO_PARCELS}"
    "OUT_ROOT=${OUT_ROOT}"
    "TP_GLOB=${row_tp_glob}"
    "THREADS=${THREADS}"
    "TEMPLATE_N4=${TEMPLATE_N4}"
    "DRIVER_PREFIX_PRIORITY=${DRIVER_PREFIX_PRIORITY}"
    "DRIVER_EXCLUDE_PREFIXES=${DRIVER_EXCLUDE_PREFIXES}"
  )
  if [[ -n ${ATLAS_TEMPLATE} ]]; then
    env_items+=("ATLAS_TEMPLATE=${ATLAS_TEMPLATE}")
  fi
  if [[ ${NO_PARCELS} -eq 0 ]]; then
    env_items+=("ATLAS_PARCELS=${ATLAS_PARCELS}")
  fi
  if [[ -n ${ATLAS_LABELS} ]]; then
    env_items+=("ATLAS_LABELS=${ATLAS_LABELS}")
  fi
  export_arg="ALL,$(IFS=,; echo "${env_items[*]}")"

  job_name="mlong_${mouse_id}"
  sbatch_cmd=(
    sbatch
    -J "${job_name}"
    -p "${PARTITION}"
    -c "${CPUS}"
    --time "${TIME_LIMIT}"
    --export "${export_arg}"
    "${SCRIPT_DIR}/run_mouse_longitudinal_pipeline.sbatch"
  )

  if [[ ${DRY_RUN} -eq 1 ]]; then
    echo "DRYRUN ${sbatch_cmd[*]}"
    printf '%s\t%s\t%s\t\tDRY_RUN\n' "${mouse_id}" "${mouse_root}" "${row_tp_glob}" >> "${SUBMIT_TSV}"
    submit_count=$((submit_count + 1))
    continue
  fi

  out="$("${sbatch_cmd[@]}")"
  job_id=$(echo "${out}" | awk '{print $NF}')
  if [[ -z ${job_id} ]]; then
    echo "WARN: could not parse job id for ${mouse_id}; sbatch output: ${out}" >&2
    printf '%s\t%s\t%s\t\tSUBMIT_PARSE_FAIL\n' "${mouse_id}" "${mouse_root}" "${row_tp_glob}" >> "${SUBMIT_TSV}"
    skip_count=$((skip_count + 1))
    continue
  fi

  echo "Submitted ${mouse_id} -> ${job_id}"
  printf '%s\t%s\t%s\t%s\tSUBMITTED\n' "${mouse_id}" "${mouse_root}" "${row_tp_glob}" "${job_id}" >> "${SUBMIT_TSV}"
  submit_count=$((submit_count + 1))
done < "${COHORT_TSV}"

echo "Done: submitted=${submit_count} skipped=${skip_count}"
echo "Submission manifest: ${SUBMIT_TSV}"
