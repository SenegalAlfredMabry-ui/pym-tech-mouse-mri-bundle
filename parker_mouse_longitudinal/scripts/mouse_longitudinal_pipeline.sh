#!/bin/bash -l
# Mouse longitudinal pipeline:
# 1) pick high-res driver per timepoint
# 2) build within-mouse unbiased template (ANTs)
# 3) warp template to Allen atlas
# 4) per-timepoint registrations with ANTs and/or @SSwarper
# 5) output transformed images, transforms, DBM Jacobians, parcel tables, and QC report

set -euo pipefail

usage() {
  cat <<'EOF'
Usage: Scripts/mouse_longitudinal_pipeline.sh [options]

Required:
  --mouse-root <dir>        Mouse folder containing timepoint dirs (e.g., MRI_data/Raw/nonhuman_mouse/1118)

Optional:
  --atlas-parcels <nii>     Allen parcel label image (integer labels)
  --no-parcels              Skip parcel warps + parcel/QC tables (alignment/Jacobians only)
  --atlas-template <nii>    Allen template image in target atlas space
                           (default: MRI_data/Raw/nonhuman_mouse/atlas/average_template_50_brain.nii.gz)
  --mouse-id <id>           Mouse ID (default: basename of --mouse-root)
  --out-root <dir>          Output root (default: MRI_data/processed/nonhuman_longitudinal)
  --methods <csv>           ants,sswarper,or both (default: ants,sswarper)
  --driver-prefix-priority <csv>
                           Preferred driver series prefixes (default: 3000,9000,4000)
  --driver-exclude-prefixes <csv>
                           Driver prefixes to avoid when alternatives exist (default: 6000,7000,8000)
  --atlas-labels <tsv/csv>  Optional parcel labels table (parcel_id + label columns)
  --timepoint-glob <glob>   Timepoint dir glob under mouse root (default: *)
  --threads <n>             Threads (default: SLURM_CPUS_PER_TASK or 4)
  --template-n4 <0|1>       Use N4 inside template build (default: 1)
  --skip-template           Reuse existing within-subject template outputs
  --skip-template2allen     Reuse existing template-to-Allen outputs
  --help

Notes:
  - Uses ANTs from: Software/ants-2.6.2/bin (same convention as ANTMAN/antpunk).
  - Loads AFNI via module (gnu13/13.2.0 + afni/25.0.11) for @SSwarper/3dNwarp tools.
EOF
}

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(cd -- "${SCRIPT_DIR}/.." && pwd)
DEFAULT_ATLAS_TEMPLATE="${REPO_ROOT}/MRI_data/Raw/nonhuman_mouse/atlas/average_template_50_brain.nii.gz"

MOUSE_ROOT=""
MOUSE_ID=""
ATLAS_TEMPLATE=""
ATLAS_PARCELS=""
ATLAS_LABELS=""
OUT_ROOT="${REPO_ROOT}/MRI_data/processed/nonhuman_longitudinal"
METHODS="ants,sswarper"
DRIVER_PREFIX_PRIORITY="3000,9000,4000"
DRIVER_EXCLUDE_PREFIXES="6000,7000,8000"
TP_GLOB="*"
THREADS="${SLURM_CPUS_PER_TASK:-4}"
TEMPLATE_N4=1
SKIP_TEMPLATE=0
SKIP_TEMPLATE2ALLEN=0
NO_PARCELS=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --mouse-root) MOUSE_ROOT=$2; shift 2 ;;
    --mouse-id) MOUSE_ID=$2; shift 2 ;;
    --atlas-template) ATLAS_TEMPLATE=$2; shift 2 ;;
    --atlas-parcels) ATLAS_PARCELS=$2; shift 2 ;;
    --no-parcels) NO_PARCELS=1; shift ;;
    --atlas-labels) ATLAS_LABELS=$2; shift 2 ;;
    --out-root) OUT_ROOT=$2; shift 2 ;;
    --methods) METHODS=$2; shift 2 ;;
    --driver-prefix-priority) DRIVER_PREFIX_PRIORITY=$2; shift 2 ;;
    --driver-exclude-prefixes) DRIVER_EXCLUDE_PREFIXES=$2; shift 2 ;;
    --timepoint-glob) TP_GLOB=$2; shift 2 ;;
    --threads) THREADS=$2; shift 2 ;;
    --template-n4) TEMPLATE_N4=$2; shift 2 ;;
    --skip-template) SKIP_TEMPLATE=1; shift ;;
    --skip-template2allen) SKIP_TEMPLATE2ALLEN=1; shift ;;
    --help|-h) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage; exit 1 ;;
  esac
done

[[ -n ${MOUSE_ROOT} ]] || { echo "ERROR: --mouse-root is required" >&2; exit 1; }
if [[ ${NO_PARCELS} -eq 0 ]]; then
  [[ -n ${ATLAS_PARCELS} ]] || { echo "ERROR: --atlas-parcels is required unless --no-parcels is set" >&2; exit 1; }
fi

if [[ -z ${ATLAS_TEMPLATE} && -f ${DEFAULT_ATLAS_TEMPLATE} ]]; then
  ATLAS_TEMPLATE=${DEFAULT_ATLAS_TEMPLATE}
  echo "INFO: using default atlas template: ${ATLAS_TEMPLATE}" >&2
fi
[[ -n ${ATLAS_TEMPLATE} ]] || {
  echo "ERROR: --atlas-template is required (default missing: ${DEFAULT_ATLAS_TEMPLATE})" >&2
  exit 1
}

MOUSE_ROOT=$(realpath "${MOUSE_ROOT}")
ATLAS_TEMPLATE=$(realpath "${ATLAS_TEMPLATE}")
if [[ -n ${ATLAS_PARCELS} ]]; then
  ATLAS_PARCELS=$(realpath "${ATLAS_PARCELS}")
fi
if [[ -n ${ATLAS_LABELS} ]]; then
  ATLAS_LABELS=$(realpath "${ATLAS_LABELS}")
fi

if [[ -z ${MOUSE_ID} ]]; then
  MOUSE_ID=$(basename "${MOUSE_ROOT}")
fi

if [[ ! -d ${MOUSE_ROOT} ]]; then
  echo "ERROR: mouse root missing: ${MOUSE_ROOT}" >&2
  exit 1
fi
if [[ ! -f ${ATLAS_TEMPLATE} ]]; then
  echo "ERROR: atlas template missing: ${ATLAS_TEMPLATE}" >&2
  exit 1
fi
if [[ ${NO_PARCELS} -eq 0 ]]; then
  if [[ ! -f ${ATLAS_PARCELS} ]]; then
    echo "ERROR: atlas parcels missing: ${ATLAS_PARCELS}" >&2
    exit 1
  fi
fi

log() { printf '[%s] %s\n' "$(date +'%F %T')" "$*"; }

# Environment setup mirroring ANTMAN/antpunk conventions
export ANTSPATH="${REPO_ROOT}/Software/ants-2.6.2/bin"
export PATH="${ANTSPATH}:${PATH}"
if command -v module >/dev/null 2>&1; then
  module load gnu13/13.2.0 afni/25.0.11 >/dev/null 2>&1 || true
fi

check_cmd() {
  command -v "$1" >/dev/null 2>&1 || { echo "ERROR: command not found: $1" >&2; exit 1; }
}
check_cmd antsRegistration
check_cmd antsApplyTransforms
check_cmd antsMultivariateTemplateConstruction2.sh
check_cmd CreateJacobianDeterminantImage
check_cmd 3dNwarpApply
check_cmd 3dNwarpFuncs
check_cmd 3dNwarpCat
check_cmd 3dcalc
check_cmd 3dTcat
check_cmd @SSwarper
check_cmd python
if [[ ${NO_PARCELS} -eq 0 ]]; then
  check_cmd LabelOverlapMeasures
fi

OUT_MOUSE="${OUT_ROOT}/${MOUSE_ID}"
INPUT_DIR="${OUT_MOUSE}/inputs"
TPL_DIR="${OUT_MOUSE}/template"
ANTS_DIR="${OUT_MOUSE}/ants"
SSW_DIR="${OUT_MOUSE}/sswarper"
TABLE_DIR="${OUT_MOUSE}/tables"
LOG_DIR="${OUT_MOUSE}/logs"
mkdir -p "${INPUT_DIR}" "${TPL_DIR}" "${ANTS_DIR}" "${SSW_DIR}" "${TABLE_DIR}" "${LOG_DIR}"

SELECTED_TSV="${INPUT_DIR}/selected_driver_scans.tsv"
log "Selecting driver scans for ${MOUSE_ID}"
python "${SCRIPT_DIR}/select_mouse_driver_scans.py" \
  --mouse-root "${MOUSE_ROOT}" \
  --timepoint-glob "${TP_GLOB}" \
  --series-prefix-priority "${DRIVER_PREFIX_PRIORITY}" \
  --exclude-prefixes "${DRIVER_EXCLUDE_PREFIXES}" \
  --output-tsv "${SELECTED_TSV}"

mapfile -t SELECTED_ROWS < <(tail -n +2 "${SELECTED_TSV}")
if [[ ${#SELECTED_ROWS[@]} -lt 2 ]]; then
  echo "ERROR: need at least 2 selected timepoints for unbiased template; found ${#SELECTED_ROWS[@]}" >&2
  exit 1
fi

TIMEPOINTS=()
DRIVERS=()
DRIVER_REF_DIR="${INPUT_DIR}/driver_reference_index"
DRIVER_REF_TSV="${INPUT_DIR}/driver_reference_manifest.tsv"
mkdir -p "${DRIVER_REF_DIR}"
printf 'timepoint\treference\tmodality_hint\tdriver_nifti\n' >"${DRIVER_REF_TSV}"

for row in "${SELECTED_ROWS[@]}"; do
  tp=$(echo "${row}" | awk -F'\t' '{print $1}')
  img=$(echo "${row}" | awk -F'\t' '{print $2}')
  modality_hint=$(echo "${row}" | awk -F'\t' '{print $16}')
  ref_name=$(echo "${row}" | awk -F'\t' '{print $17}')
  if [[ -z ${ref_name} ]]; then
    ref_name="${MOUSE_ID}_${tp}_driver"
  fi
  safe_ref=$(echo "${ref_name}" | sed 's/[^A-Za-z0-9._-]/_/g')
  img_use="${img}"
  n4=$(3dinfo -n4 "${img}" 2>/dev/null | awk '{print $NF}' || echo 1)
  if [[ ! ${n4} =~ ^[0-9]+$ ]]; then
    n4=1
  fi
  if [[ ${n4} -gt 1 ]]; then
    img_use="${INPUT_DIR}/${tp}_driver_3d.nii.gz"
    log "Timepoint ${tp}: converting 4D driver (${n4} vols) -> 3D first volume"
    3dTcat -overwrite -prefix "${img_use}" "${img}[0]" >/dev/null
  fi
  TIMEPOINTS+=("${tp}")
  DRIVERS+=("${img_use}")
  ln -sfn "${img_use}" "${DRIVER_REF_DIR}/${safe_ref}.nii.gz"
  printf '%s\t%s\t%s\t%s\n' "${tp}" "${safe_ref}" "${modality_hint}" "${img_use}" >>"${DRIVER_REF_TSV}"
done

log "Selected ${#TIMEPOINTS[@]} timepoints: ${TIMEPOINTS[*]}"
log "Driver reference index: ${DRIVER_REF_DIR}"

TEMPLATE_PREFIX="${TPL_DIR}/${MOUSE_ID}_within_"
RAW_TEMPLATE="${TEMPLATE_PREFIX}template0.nii.gz"
WITHIN_TEMPLATE="${TPL_DIR}/within_subject_template.nii.gz"

if [[ ${SKIP_TEMPLATE} -eq 0 || ! -f ${WITHIN_TEMPLATE} ]]; then
  log "Building unbiased within-subject template (ANTs)"
  antsMultivariateTemplateConstruction2.sh \
    -d 3 \
    -o "${TEMPLATE_PREFIX}" \
    -i 4 \
    -c 2 \
    -j "${THREADS}" \
    -n "${TEMPLATE_N4}" \
    -k 1 \
    -w 1 \
    -g 0.2 \
    -f 6x4x2x1 \
    -s 3x2x1x0vox \
    -q 120x80x40x20 \
    -m CC[4] \
    -t SyN[0.08,3,0] \
    "${DRIVERS[@]}" \
    >"${LOG_DIR}/template_build_${MOUSE_ID}.log" 2>&1

  if [[ ! -f ${RAW_TEMPLATE} ]]; then
    echo "ERROR: expected template not found: ${RAW_TEMPLATE}" >&2
    echo "See log: ${LOG_DIR}/template_build_${MOUSE_ID}.log" >&2
    exit 1
  fi
  cp -f "${RAW_TEMPLATE}" "${WITHIN_TEMPLATE}"
else
  log "Reusing existing within-subject template: ${WITHIN_TEMPLATE}"
fi

TPL2ALLEN_PREFIX="${TPL_DIR}/template_to_allen_"
TPL_IN_ALLEN="${TPL_DIR}/within_template_in_allen.nii.gz"
ALLEN_IN_TPL="${TPL_DIR}/allen_template_in_within.nii.gz"

if [[ ${SKIP_TEMPLATE2ALLEN} -eq 0 || ! -f ${TPL2ALLEN_PREFIX}1Warp.nii.gz || ! -f ${TPL2ALLEN_PREFIX}1InverseWarp.nii.gz ]]; then
  log "Registering within-subject template to Allen template"
  antsRegistration \
    --dimensionality 3 \
    --float 1 \
    --output "[${TPL2ALLEN_PREFIX},${TPL_IN_ALLEN},${ALLEN_IN_TPL}]" \
    --interpolation Linear \
    --winsorize-image-intensities "[0.005,0.995]" \
    --use-histogram-matching 1 \
    --initial-moving-transform "[${ATLAS_TEMPLATE},${WITHIN_TEMPLATE},1]" \
    --transform Rigid[0.08] \
    --metric MI[${ATLAS_TEMPLATE},${WITHIN_TEMPLATE},1,32,Regular,0.25] \
    --convergence "[1000x500x250x0,1e-6,10]" \
    --shrink-factors 8x4x2x1 \
    --smoothing-sigmas 3x2x1x0vox \
    --transform Affine[0.08] \
    --metric MI[${ATLAS_TEMPLATE},${WITHIN_TEMPLATE},1,32,Regular,0.25] \
    --convergence "[1000x500x250x0,1e-6,10]" \
    --shrink-factors 8x4x2x1 \
    --smoothing-sigmas 3x2x1x0vox \
    --transform SyN[0.08,3,0] \
    --metric CC[${ATLAS_TEMPLATE},${WITHIN_TEMPLATE},1,4] \
    --convergence "[120x80x40x20,1e-7,10]" \
    --shrink-factors 8x4x2x1 \
    --smoothing-sigmas 3x2x1x0vox \
    >"${LOG_DIR}/template_to_allen_${MOUSE_ID}.log" 2>&1
fi

if [[ ! -f ${TPL2ALLEN_PREFIX}1Warp.nii.gz || ! -f ${TPL2ALLEN_PREFIX}1InverseWarp.nii.gz || ! -f ${TPL2ALLEN_PREFIX}0GenericAffine.mat ]]; then
  echo "ERROR: template-to-Allen transforms missing" >&2
  exit 1
fi

TPL_PARCELS="${TPL_DIR}/allen_parcels_in_within_template.nii.gz"
if [[ ${NO_PARCELS} -eq 0 ]]; then
  log "Projecting Allen parcels into within-subject template space"
  antsApplyTransforms -d 3 \
    -i "${ATLAS_PARCELS}" \
    -r "${WITHIN_TEMPLATE}" \
    -n NearestNeighbor \
    -t "[${TPL2ALLEN_PREFIX}0GenericAffine.mat,1]" \
    -t "${TPL2ALLEN_PREFIX}1InverseWarp.nii.gz" \
    -o "${TPL_PARCELS}" \
    >"${LOG_DIR}/parcels_to_template_${MOUSE_ID}.log" 2>&1
else
  log "Skipping parcel projection (--no-parcels)"
fi

MANIFEST="${TABLE_DIR}/run_manifest.tsv"
TRANSFORM_MANIFEST="${TABLE_DIR}/transforms_manifest.tsv"
: >"${MANIFEST}"
: >"${TRANSFORM_MANIFEST}"
printf 'method\ttimepoint\timage_allen\tjacobian\troundtrip_parcels\n' >"${MANIFEST}"
printf 'method\ttimepoint\tnative_image\timage_template\tforward_transform_1\tforward_transform_2\n' >"${TRANSFORM_MANIFEST}"

has_method() {
  local needle=$1
  IFS=',' read -ra _m <<<"${METHODS}"
  for m in "${_m[@]}"; do
    if [[ $(echo "$m" | xargs) == "${needle}" ]]; then
      return 0
    fi
  done
  return 1
}

run_ants_tp() {
  local tp=$1
  local native_img=$2
  local out_tp="${ANTS_DIR}/${tp}"
  mkdir -p "${out_tp}"

  local pref="${out_tp}/${tp}_to_template_"
  local img_tpl="${out_tp}/${tp}_in_template.nii.gz"
  local img_allen="${out_tp}/${tp}_in_allen.nii.gz"
  local jac="${out_tp}/${tp}_jacobian_template.nii.gz"
  local native_parc="${out_tp}/${tp}_parcels_native_from_template.nii.gz"
  local roundtrip="${out_tp}/${tp}_parcels_roundtrip_template.nii.gz"

  log "[ANTs] ${tp}: registering native -> within template"
  antsRegistration \
    --dimensionality 3 \
    --float 1 \
    --output "[${pref},${img_tpl}]" \
    --interpolation Linear \
    --winsorize-image-intensities "[0.005,0.995]" \
    --use-histogram-matching 1 \
    --initial-moving-transform "[${WITHIN_TEMPLATE},${native_img},1]" \
    --transform Rigid[0.08] \
    --metric MI[${WITHIN_TEMPLATE},${native_img},1,32,Regular,0.25] \
    --convergence "[1000x500x250x0,1e-6,10]" \
    --shrink-factors 8x4x2x1 \
    --smoothing-sigmas 3x2x1x0vox \
    --transform Affine[0.08] \
    --metric MI[${WITHIN_TEMPLATE},${native_img},1,32,Regular,0.25] \
    --convergence "[1000x500x250x0,1e-6,10]" \
    --shrink-factors 8x4x2x1 \
    --smoothing-sigmas 3x2x1x0vox \
    --transform SyN[0.08,3,0] \
    --metric CC[${WITHIN_TEMPLATE},${native_img},1,4] \
    --convergence "[120x80x40x20,1e-7,10]" \
    --shrink-factors 8x4x2x1 \
    --smoothing-sigmas 3x2x1x0vox \
    >"${out_tp}/${tp}_ants_registration.log" 2>&1

  CreateJacobianDeterminantImage 3 "${pref}1Warp.nii.gz" "${jac}" 1 0

  antsApplyTransforms -d 3 \
    -i "${native_img}" \
    -r "${ATLAS_TEMPLATE}" \
    -o "${img_allen}" \
    -n Linear \
    -t "${TPL2ALLEN_PREFIX}1Warp.nii.gz" \
    -t "${TPL2ALLEN_PREFIX}0GenericAffine.mat" \
    -t "${pref}1Warp.nii.gz" \
    -t "${pref}0GenericAffine.mat"

  if [[ ${NO_PARCELS} -eq 0 ]]; then
    # Roundtrip parcels for Dice/landmark QC.
    antsApplyTransforms -d 3 \
      -i "${TPL_PARCELS}" \
      -r "${native_img}" \
      -o "${native_parc}" \
      -n NearestNeighbor \
      -t "[${pref}0GenericAffine.mat,1]" \
      -t "${pref}1InverseWarp.nii.gz"

    antsApplyTransforms -d 3 \
      -i "${native_parc}" \
      -r "${WITHIN_TEMPLATE}" \
      -o "${roundtrip}" \
      -n NearestNeighbor \
      -t "${pref}1Warp.nii.gz" \
      -t "${pref}0GenericAffine.mat"

    printf 'ants\t%s\t%s\t%s\t%s\n' "${tp}" "${img_allen}" "${jac}" "${roundtrip}" >>"${MANIFEST}"
  else
    printf 'ants\t%s\t%s\t%s\t\n' "${tp}" "${img_allen}" "${jac}" >>"${MANIFEST}"
  fi
  printf 'ants\t%s\t%s\t%s\t%s\t%s\n' "${tp}" "${native_img}" "${img_tpl}" "${pref}1Warp.nii.gz" "${pref}0GenericAffine.mat" >>"${TRANSFORM_MANIFEST}"
}

pick_ss_output() {
  local base=$1
  if [[ -f ${base}.nii.gz ]]; then
    echo "${base}.nii.gz"
  elif [[ -f ${base}.nii ]]; then
    echo "${base}.nii"
  else
    echo ""
  fi
}

run_sswarper_tp() {
  local tp=$1
  local native_img=$2
  local out_tp="${SSW_DIR}/${tp}"
  mkdir -p "${out_tp}"

  log "[SSwarper] ${tp}: registering native -> within template"
  @SSwarper \
    -input "${native_img}" \
    -base "${WITHIN_TEMPLATE}" \
    -subid "${tp}" \
    -odir "${out_tp}" \
    -warpscale 0.8 \
    -minp 9 \
    >"${out_tp}/${tp}_sswarper.log" 2>&1

  local img_tpl
  img_tpl=$(pick_ss_output "${out_tp}/anatQQ.${tp}")
  if [[ -z ${img_tpl} ]]; then
    echo "ERROR: sswarper output anatQQ missing for ${tp}" >&2
    exit 1
  fi

  local warp_field="${out_tp}/anatQQ.${tp}_WARP.nii"
  local warp_aff="${out_tp}/anatQQ.${tp}.aff12.1D"
  if [[ ! -f ${warp_field} && -f ${warp_field}.gz ]]; then
    warp_field="${warp_field}.gz"
  fi
  if [[ ! -f ${warp_field} || ! -f ${warp_aff} ]]; then
    echo "ERROR: sswarper transforms missing for ${tp}" >&2
    exit 1
  fi

  local standardized_tpl="${out_tp}/${tp}_in_template.nii.gz"
  cp -f "${img_tpl}" "${standardized_tpl}"

  local forward_comp="${out_tp}/${tp}_to_template_composite_WARP.nii.gz"
  3dNwarpCat -overwrite -prefix "${forward_comp}" "${warp_field}" "${warp_aff}" >/dev/null

  local bulk="${out_tp}/${tp}_bulk_template.nii.gz"
  local jac="${out_tp}/${tp}_jacobian_template.nii.gz"
  3dNwarpFuncs -nwarp "${forward_comp}" -bulk -prefix "${bulk}" >/dev/null
  3dcalc -a "${bulk}" -expr 'a+1' -prefix "${jac}" >/dev/null

  local img_allen="${out_tp}/${tp}_in_allen.nii.gz"
  antsApplyTransforms -d 3 \
    -i "${standardized_tpl}" \
    -r "${ATLAS_TEMPLATE}" \
    -o "${img_allen}" \
    -n Linear \
    -t "${TPL2ALLEN_PREFIX}1Warp.nii.gz" \
    -t "${TPL2ALLEN_PREFIX}0GenericAffine.mat"

  if [[ ${NO_PARCELS} -eq 0 ]]; then
    local native_parc="${out_tp}/${tp}_parcels_native_from_template.nii.gz"
    local roundtrip="${out_tp}/${tp}_parcels_roundtrip_template.nii.gz"
    3dNwarpApply -overwrite \
      -master "${native_img}" \
      -nwarp "INV(${forward_comp})" \
      -source "${TPL_PARCELS}" \
      -ainterp NN \
      -prefix "${native_parc}" >/dev/null

    3dNwarpApply -overwrite \
      -master "${WITHIN_TEMPLATE}" \
      -nwarp "${forward_comp}" \
      -source "${native_parc}" \
      -ainterp NN \
      -prefix "${roundtrip}" >/dev/null

    printf 'sswarper\t%s\t%s\t%s\t%s\n' "${tp}" "${img_allen}" "${jac}" "${roundtrip}" >>"${MANIFEST}"
  else
    printf 'sswarper\t%s\t%s\t%s\t\n' "${tp}" "${img_allen}" "${jac}" >>"${MANIFEST}"
  fi
  printf 'sswarper\t%s\t%s\t%s\t%s\t%s\n' "${tp}" "${native_img}" "${standardized_tpl}" "${forward_comp}" "${warp_aff}" >>"${TRANSFORM_MANIFEST}"
}

for i in "${!TIMEPOINTS[@]}"; do
  tp="${TIMEPOINTS[$i]}"
  native_img="${DRIVERS[$i]}"

  if has_method ants; then
    run_ants_tp "${tp}" "${native_img}"
  fi
  if has_method sswarper; then
    run_sswarper_tp "${tp}" "${native_img}"
  fi
done

if [[ ${NO_PARCELS} -eq 0 ]]; then
  log "Aggregating parcel/Jacobian/QC tables"
  METRICS_CMD=(
    python "${SCRIPT_DIR}/mouse_longitudinal_metrics.py"
    --manifest "${MANIFEST}"
    --allen-parcels "${ATLAS_PARCELS}"
    --template-parcels "${TPL_PARCELS}"
    --out-dir "${TABLE_DIR}"
  )
  if [[ -n ${ATLAS_LABELS} && -f ${ATLAS_LABELS} ]]; then
    METRICS_CMD+=(--atlas-labels "${ATLAS_LABELS}")
  fi
  "${METRICS_CMD[@]}"
else
  log "Skipping parcel metrics/QC tables (--no-parcels)"
  printf 'Parcel-level metrics were skipped with --no-parcels on %s\n' "$(date +'%F %T')" >"${TABLE_DIR}/PARCELS_SKIPPED.txt"
fi

log "Done. Outputs: ${OUT_MOUSE}"
