#!/bin/bash -l
#SBATCH -J mouse_long
#SBATCH -p cloud
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --time=0
#SBATCH --mem=96G
#SBATCH -o sbatch-mouse-long-%j.out
#SBATCH -e sbatch-mouse-long-%j.err

set -euo pipefail

# Required env vars:
#   MOUSE_ROOT
# Optional:
#   ATLAS_PARCELS      (required unless NO_PARCELS=1)
#   ATLAS_TEMPLATE     (default: MRI_data/Raw/nonhuman_mouse/atlas/average_template_50_brain.nii.gz)
#   MOUSE_ID
#   OUT_ROOT
#   METHODS           (default: ants,sswarper)
#   DRIVER_PREFIX_PRIORITY   (default: 3000,9000,4000)
#   DRIVER_EXCLUDE_PREFIXES  (default: 6000,7000,8000)
#   ATLAS_LABELS
#   TP_GLOB           (default: *)
#   THREADS           (default: SLURM_CPUS_PER_TASK)
#   TEMPLATE_N4       (0/1; default: 1)
#   SKIP_TEMPLATE     (0/1)
#   SKIP_TEMPLATE2ALLEN (0/1)
#   NO_PARCELS        (0/1; if 1, run alignment/Jacobians only)

REPO_ROOT=${SNHEART_ROOT:-/path/to/SNHEART}
DEFAULT_ATLAS_TEMPLATE=${REPO_ROOT}/MRI_data/Raw/nonhuman_mouse/atlas/average_template_50_brain.nii.gz
cd "${REPO_ROOT}"

: "${MOUSE_ROOT:?Set MOUSE_ROOT}"
ATLAS_TEMPLATE=${ATLAS_TEMPLATE:-${DEFAULT_ATLAS_TEMPLATE}}
NO_PARCELS=${NO_PARCELS:-0}

if [[ ${NO_PARCELS} -ne 1 ]]; then
  : "${ATLAS_PARCELS:?Set ATLAS_PARCELS (or set NO_PARCELS=1)}"
fi

if [[ ! -f ${ATLAS_TEMPLATE} ]]; then
  echo "ERROR: atlas template missing: ${ATLAS_TEMPLATE}" >&2
  exit 1
fi

METHODS=${METHODS:-ants,sswarper}
DRIVER_PREFIX_PRIORITY=${DRIVER_PREFIX_PRIORITY:-3000,9000,4000}
DRIVER_EXCLUDE_PREFIXES=${DRIVER_EXCLUDE_PREFIXES:-6000,7000,8000}
OUT_ROOT=${OUT_ROOT:-${REPO_ROOT}/MRI_data/processed/nonhuman_longitudinal}
TP_GLOB=${TP_GLOB:-*}
THREADS=${THREADS:-${SLURM_CPUS_PER_TASK:-8}}
TEMPLATE_N4=${TEMPLATE_N4:-1}
SKIP_TEMPLATE=${SKIP_TEMPLATE:-0}
SKIP_TEMPLATE2ALLEN=${SKIP_TEMPLATE2ALLEN:-0}

cmd=(
  "${REPO_ROOT}/Scripts/mouse_longitudinal_pipeline.sh"
  --mouse-root "${MOUSE_ROOT}"
  --atlas-template "${ATLAS_TEMPLATE}"
  --out-root "${OUT_ROOT}"
  --methods "${METHODS}"
  --driver-prefix-priority "${DRIVER_PREFIX_PRIORITY}"
  --driver-exclude-prefixes "${DRIVER_EXCLUDE_PREFIXES}"
  --timepoint-glob "${TP_GLOB}"
  --threads "${THREADS}"
  --template-n4 "${TEMPLATE_N4}"
)

if [[ ${NO_PARCELS} -eq 1 ]]; then
  cmd+=(--no-parcels)
else
  cmd+=(--atlas-parcels "${ATLAS_PARCELS}")
fi

if [[ -n ${MOUSE_ID:-} ]]; then
  cmd+=(--mouse-id "${MOUSE_ID}")
fi
if [[ -n ${ATLAS_LABELS:-} ]]; then
  cmd+=(--atlas-labels "${ATLAS_LABELS}")
fi
if [[ ${SKIP_TEMPLATE} -eq 1 ]]; then
  cmd+=(--skip-template)
fi
if [[ ${SKIP_TEMPLATE2ALLEN} -eq 1 ]]; then
  cmd+=(--skip-template2allen)
fi

printf 'Running command:\n%s\n' "${cmd[*]}"
"${cmd[@]}"

printf 'Completed at %s\n' "$(date +'%F %T')"
