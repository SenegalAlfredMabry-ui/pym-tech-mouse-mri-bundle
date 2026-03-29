#!/bin/bash -l
# Integrated multimodal longitudinal alignment:
# Reuses the existing anatomical longitudinal transform chain and adds:
# - one T2-weighted scan per timepoint
# - all fMRI 6000/7000 scans per timepoint
#
# Transform chain used for each modality image:
# modality -> structural(native, same timepoint) -> within-template -> Allen

set -euo pipefail

usage() {
  cat <<'EOF'
Usage: Scripts/mouse_longitudinal_fmri_t2tse.sh [options]

Required:
  --mouse-root <dir>      Mouse folder with timepoint dirs (e.g., .../NIFTI_longitudinal/100M1)
  --anat-run-dir <dir>    Existing anatomical longitudinal output root for this mouse
                          (contains tables/transforms_manifest.tsv and template/)

Optional:
  --mouse-id <id>         Output tag (default: basename of --mouse-root)
  --out-root <dir>        Output root
                          (default: <anat-run-dir>/modalities_fmri_t2tse)
  --timepoint-glob <glob> Timepoint dir glob under --mouse-root (default: *)
  --threads <n>           ITK threads for registrations (default: SLURM_CPUS_PER_TASK or 4)
  --modality-nonlinear <0|1>
                          Add SyN stage for modality->structural step (default: 0)
  --help

Notes:
  - This script does NOT rebuild the longitudinal template.
  - It consumes anatomical transforms from the ANTs rows in:
      <anat-run-dir>/tables/transforms_manifest.tsv
  - T2 selection per timepoint:
      prefer *_T2tse_40001_*.nii.gz, then *_T2tse_*.nii.gz,
      then T2_RARE/T2RARE fallbacks (for alternate naming conventions).
  - fMRI selection per timepoint: *_fMRI_6000_*.nii.gz and *_fMRI_7000_*.nii.gz.
  - 4D fMRI handling: register a 3D mean reference, then apply transforms to full 4D data.
    For low-memory safety, 4D outputs are warped volume-by-volume and re-stacked.
EOF
}

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(cd -- "${SCRIPT_DIR}/.." && pwd)

MOUSE_ROOT=""
ANAT_RUN_DIR=""
MOUSE_ID=""
OUT_ROOT=""
TP_GLOB="*"
THREADS="${SLURM_CPUS_PER_TASK:-4}"
MODALITY_NONLINEAR=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --mouse-root) MOUSE_ROOT=$2; shift 2 ;;
    --anat-run-dir) ANAT_RUN_DIR=$2; shift 2 ;;
    --mouse-id) MOUSE_ID=$2; shift 2 ;;
    --out-root) OUT_ROOT=$2; shift 2 ;;
    --timepoint-glob) TP_GLOB=$2; shift 2 ;;
    --threads) THREADS=$2; shift 2 ;;
    --modality-nonlinear) MODALITY_NONLINEAR=$2; shift 2 ;;
    --help|-h) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage; exit 1 ;;
  esac
done

[[ -n ${MOUSE_ROOT} ]] || { echo "ERROR: --mouse-root is required" >&2; exit 1; }
[[ -n ${ANAT_RUN_DIR} ]] || { echo "ERROR: --anat-run-dir is required" >&2; exit 1; }
MOUSE_ROOT=$(realpath "${MOUSE_ROOT}")
ANAT_RUN_DIR=$(realpath "${ANAT_RUN_DIR}")
[[ -d ${MOUSE_ROOT} ]] || { echo "ERROR: mouse root missing: ${MOUSE_ROOT}" >&2; exit 1; }
[[ -d ${ANAT_RUN_DIR} ]] || { echo "ERROR: anatomical run dir missing: ${ANAT_RUN_DIR}" >&2; exit 1; }

if [[ -z ${MOUSE_ID} ]]; then
  MOUSE_ID=$(basename "${MOUSE_ROOT}")
fi
if [[ -z ${OUT_ROOT} ]]; then
  OUT_ROOT="${ANAT_RUN_DIR}/modalities_fmri_t2tse"
fi
OUT_ROOT=$(realpath -m "${OUT_ROOT}")

if [[ ! ${THREADS} =~ ^[0-9]+$ || ${THREADS} -lt 1 ]]; then
  echo "ERROR: --threads must be a positive integer" >&2
  exit 1
fi
if [[ ! ${MODALITY_NONLINEAR} =~ ^[01]$ ]]; then
  echo "ERROR: --modality-nonlinear must be 0 or 1" >&2
  exit 1
fi

log() { printf '[%s] %s\n' "$(date +'%F %T')" "$*"; }

# Environment setup consistent with the anatomical script.
export ANTSPATH="${REPO_ROOT}/Software/ants-2.6.2/bin"
export PATH="${ANTSPATH}:${PATH}"
export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS="${THREADS}"
if command -v module >/dev/null 2>&1; then
  module load gnu13/13.2.0 afni/25.0.11 >/dev/null 2>&1 || true
fi

check_cmd() {
  command -v "$1" >/dev/null 2>&1 || { echo "ERROR: command not found: $1" >&2; exit 1; }
}
check_cmd antsRegistration
check_cmd antsApplyTransforms
check_cmd 3dinfo
check_cmd 3dTstat
check_cmd 3dTcat
check_cmd awk
check_cmd find
check_cmd sort
check_cmd python

TRANSFORM_MANIFEST_IN="${ANAT_RUN_DIR}/tables/transforms_manifest.tsv"
WITHIN_TEMPLATE="${ANAT_RUN_DIR}/template/within_subject_template.nii.gz"
TPL2ALLEN_WARP="${ANAT_RUN_DIR}/template/template_to_allen_1Warp.nii.gz"
TPL2ALLEN_AFFINE="${ANAT_RUN_DIR}/template/template_to_allen_0GenericAffine.mat"
ATLAS_TEMPLATE="${ANAT_RUN_DIR}/template/within_template_in_allen.nii.gz"

[[ -f ${TRANSFORM_MANIFEST_IN} ]] || { echo "ERROR: missing ${TRANSFORM_MANIFEST_IN}" >&2; exit 1; }
[[ -f ${WITHIN_TEMPLATE} ]] || { echo "ERROR: missing ${WITHIN_TEMPLATE}" >&2; exit 1; }
[[ -f ${TPL2ALLEN_WARP} ]] || { echo "ERROR: missing ${TPL2ALLEN_WARP}" >&2; exit 1; }
[[ -f ${TPL2ALLEN_AFFINE} ]] || { echo "ERROR: missing ${TPL2ALLEN_AFFINE}" >&2; exit 1; }
[[ -f ${ATLAS_TEMPLATE} ]] || { echo "ERROR: missing ${ATLAS_TEMPLATE}" >&2; exit 1; }

hash_driver_file() {
  local file=$1
  if command -v sha256sum >/dev/null 2>&1; then
    sha256sum "${file}" | awk '{print $1}'
  elif command -v md5sum >/dev/null 2>&1; then
    md5sum "${file}" | awk '{print $1}'
  else
    stat -c '%Y%z%s' "${file}"
  fi
}

prepare_driver_volumes() {
  local driver=$1
  local hash
  hash=$(hash_driver_file "${driver}")
  local base
  base=$(basename "${driver}" .nii.gz)
  local prep_dir="${FMRI_DRIVER_PREP_DIR}/${hash}_${base}"
  local prep_tar="${prep_dir}.tar.gz"
  local count_file="${prep_dir}/volume_count.txt"

  local nvol
  nvol=$(3dinfo -nv "${driver}" 2>/dev/null | awk '{print $1}' || echo 0)
  if [[ ! ${nvol} =~ ^[0-9]+$ || ${nvol} -lt 1 ]]; then
    echo "ERROR: unable to determine volumes for ${driver}" >&2
    exit 1
  fi

  if [[ -f ${prep_tar} && ! -d ${prep_dir} ]]; then
    mkdir -p "${prep_dir}"
    tar -xzf "${prep_tar}" -C "${prep_dir}"
  fi

  if [[ -d ${prep_dir} && -f ${count_file} ]]; then
    local stored
    stored=$(<"${count_file}")
    if [[ "${stored}" -eq "${nvol}" ]]; then
      echo "${prep_dir}"
      return 0
    fi
    rm -rf "${prep_dir}" "${prep_tar}"
  fi

  mkdir -p "${prep_dir}"
  for ((i=0; i<nvol; i++)); do
    printf -v idx "%06d" "${i}"
    3dTcat -overwrite -prefix "${prep_dir}/vol_${idx}.nii.gz" "${driver}[${i}]" >/dev/null
  done

  printf '%d\n' "${nvol}" >"${count_file}"
  mkdir -p "$(dirname "${prep_tar}")"
  tar -czf "${prep_tar}" -C "$(dirname "${prep_dir}")" "$(basename "${prep_dir}")" >/dev/null
  echo "${prep_dir}"
}

OUT_MOUSE="${OUT_ROOT}/${MOUSE_ID}"
INPUT_DIR="${OUT_MOUSE}/inputs"
ALIGN_DIR="${OUT_MOUSE}/aligned"
TABLE_DIR="${OUT_MOUSE}/tables"
LOG_DIR="${OUT_MOUSE}/logs"
mkdir -p "${INPUT_DIR}" "${ALIGN_DIR}" "${TABLE_DIR}" "${LOG_DIR}"
FMRI_DRIVER_PREP_DIR="${FMRI_DRIVER_PREP_DIR:-${OUT_MOUSE}/driver_preps}"
mkdir -p "${FMRI_DRIVER_PREP_DIR}"
FMRI_TCAT_CHUNK_SIZE="${FMRI_TCAT_CHUNK_SIZE:-12}"

SELECTED_TSV="${INPUT_DIR}/selected_t2_and_fmri.tsv"
MANIFEST_TSV="${TABLE_DIR}/alignment_manifest.tsv"
CHAIN_TSV="${TABLE_DIR}/transform_chain_manifest.tsv"

printf 'timepoint\tstructural_native\tt2_selected\tn_fmri\tfmri_images\n' >"${SELECTED_TSV}"
printf 'timepoint\trole\tlabel\tsource_image\tsource_ndim\tref3d_image\tstructural_native\taligned_struct\taligned_template\taligned_allen\tregistration_log\n' >"${MANIFEST_TSV}"
printf 'timepoint\trole\tlabel\tmod_to_struct_warp\tmod_to_struct_inverse_warp\tmod_to_struct_affine\tstruct_to_template_warp\tstruct_to_template_affine\ttemplate_to_allen_warp\ttemplate_to_allen_affine\n' >"${CHAIN_TSV}"

declare -A STRUCT_NATIVE
declare -A STRUCT_TO_TEMPLATE_WARP
declare -A STRUCT_TO_TEMPLATE_AFFINE
declare -A STRUCT_TO_TEMPLATE_INVWARP

while IFS=$'\t' read -r method tp native_image image_template fwd1 fwd2 _rest; do
  [[ ${method} == "method" ]] && continue
  [[ ${method} == "ants" ]] || continue
  [[ -n ${tp} ]] || continue
  [[ -n ${native_image} && -n ${fwd1} && -n ${fwd2} ]] || continue

  STRUCT_NATIVE["${tp}"]="${native_image}"
  STRUCT_TO_TEMPLATE_WARP["${tp}"]="${fwd1}"
  STRUCT_TO_TEMPLATE_AFFINE["${tp}"]="${fwd2}"

  inv_candidate=""
  if [[ ${fwd1} == *"1Warp.nii.gz" ]]; then
    inv_candidate="${fwd1%1Warp.nii.gz}1InverseWarp.nii.gz"
  fi
  STRUCT_TO_TEMPLATE_INVWARP["${tp}"]="${inv_candidate}"
done < "${TRANSFORM_MANIFEST_IN}"

if [[ ${#STRUCT_NATIVE[@]} -eq 0 ]]; then
  echo "ERROR: no ANTs rows found in ${TRANSFORM_MANIFEST_IN}" >&2
  exit 1
fi

pick_t2_driver() {
  local tp_dir=$1
  local preferred

  # Primary preference: standard T2tse naming.
  preferred=$(find "${tp_dir}" -maxdepth 1 -type f -iname '*_T2tse_40001_*.nii.gz' | sort | head -n 1 || true)
  if [[ -n ${preferred} ]]; then
    echo "${preferred}"
    return 0
  fi

  preferred=$(find "${tp_dir}" -maxdepth 1 -type f -iname '*_T2tse_*.nii.gz' | sort | head -n 1 || true)
  if [[ -n ${preferred} ]]; then
    echo "${preferred}"
    return 0
  fi

  # Fallback naming seen in converted mouse datasets.
  preferred=$(find "${tp_dir}" -maxdepth 1 -type f -iname '*_T2_RARE_90001_*.nii.gz' | sort | head -n 1 || true)
  if [[ -n ${preferred} ]]; then
    echo "${preferred}"
    return 0
  fi

  preferred=$(find "${tp_dir}" -maxdepth 1 -type f -iname '*_T2_RARE_*.nii.gz' | sort | head -n 1 || true)
  if [[ -n ${preferred} ]]; then
    echo "${preferred}"
    return 0
  fi

  preferred=$(find "${tp_dir}" -maxdepth 1 -type f -iname '*_T2RARE_*.nii.gz' | sort | head -n 1 || true)
  if [[ -n ${preferred} ]]; then
    echo "${preferred}"
    return 0
  fi

  # Final loose fallback for unexpected case/style variants.
  find "${tp_dir}" -maxdepth 1 -type f -iname '*t2*rare*.nii.gz' | sort | head -n 1 || true
}

prepare_ref3d() {
  local src=$1
  local out_ref=$2
  local n4
  n4=$(3dinfo -n4 "${src}" 2>/dev/null | awk '{print $NF}' || echo 1)
  if [[ ! ${n4} =~ ^[0-9]+$ ]]; then
    n4=1
  fi
  if [[ ${n4} -gt 1 ]]; then
    3dTstat -overwrite -mean -prefix "${out_ref}" "${src}" >/dev/null
  else
    cp -f "${src}" "${out_ref}"
  fi
  echo "${n4}"
}

register_modality_to_struct() {
  local struct_img=$1
  local moving_ref3d=$2
  local out_prefix=$3
  local out_ref_in_struct=$4
  local reg_log=$5

  if [[ ${MODALITY_NONLINEAR} -eq 1 ]]; then
    antsRegistration \
      --dimensionality 3 \
      --float 1 \
      --output "[${out_prefix},${out_ref_in_struct}]" \
      --interpolation Linear \
      --winsorize-image-intensities "[0.005,0.995]" \
      --use-histogram-matching 1 \
      --initial-moving-transform "[${struct_img},${moving_ref3d},1]" \
      --transform Rigid[0.08] \
      --metric MI[${struct_img},${moving_ref3d},1,32,Regular,0.25] \
      --convergence "[1000x500x250x0,1e-6,10]" \
      --shrink-factors 8x4x2x1 \
      --smoothing-sigmas 3x2x1x0vox \
      --transform Affine[0.08] \
      --metric MI[${struct_img},${moving_ref3d},1,32,Regular,0.25] \
      --convergence "[1000x500x250x0,1e-6,10]" \
      --shrink-factors 8x4x2x1 \
      --smoothing-sigmas 3x2x1x0vox \
      --transform SyN[0.08,3,0] \
      --metric CC[${struct_img},${moving_ref3d},1,4] \
      --convergence "[120x80x40x20,1e-7,10]" \
      --shrink-factors 8x4x2x1 \
      --smoothing-sigmas 3x2x1x0vox \
      >"${reg_log}" 2>&1
  else
    antsRegistration \
      --dimensionality 3 \
      --float 1 \
      --output "[${out_prefix},${out_ref_in_struct}]" \
      --interpolation Linear \
      --winsorize-image-intensities "[0.005,0.995]" \
      --use-histogram-matching 1 \
      --initial-moving-transform "[${struct_img},${moving_ref3d},1]" \
      --transform Rigid[0.08] \
      --metric MI[${struct_img},${moving_ref3d},1,32,Regular,0.25] \
      --convergence "[1000x500x250x0,1e-6,10]" \
      --shrink-factors 8x4x2x1 \
      --smoothing-sigmas 3x2x1x0vox \
      --transform Affine[0.08] \
      --metric MI[${struct_img},${moving_ref3d},1,32,Regular,0.25] \
      --convergence "[1000x500x250x0,1e-6,10]" \
      --shrink-factors 8x4x2x1 \
      --smoothing-sigmas 3x2x1x0vox \
      >"${reg_log}" 2>&1
  fi
}

apply_chain() {
  local src=$1
  local src_ndim=$2
  local struct_native=$3
  local out_in_struct=$4
  local out_in_template=$5
  local out_in_allen=$6
  local mod_warp=$7
  local mod_affine=$8
  local struct_warp=$9
  local struct_affine=${10}

  apply_3d_once() {
    local in_img=$1
    local ref_img=$2
    local out_img=$3
    shift 3
    local xforms=("$@")
    antsApplyTransforms -d 3 -i "${in_img}" -r "${ref_img}" -o "${out_img}" -n Linear "${xforms[@]}"
  }

  apply_4d_lowmem() {
    local in_4d=$1
    local ref_img=$2
    local out_4d=$3
    shift 3
    local xforms=("$@")

    if [[ -f ${out_4d} ]]; then
      log "Reusing existing 4D output: $(basename "${out_4d}")"
      return 0
    fi

    local nvol
    nvol=$(3dinfo -nv "${in_4d}" 2>/dev/null | awk '{print $1}' || echo 0)
    if [[ ! ${nvol} =~ ^[0-9]+$ || ${nvol} -lt 1 ]]; then
      echo "ERROR: unable to determine number of volumes for ${in_4d}" >&2
      exit 1
    fi

    local tmp_dir
    local tmp_candidates=()
    shopt -s nullglob
    tmp_candidates=("${out_4d}.tmp."*)
    shopt -u nullglob
    if (( ${#tmp_candidates[@]} > 1 )); then
      echo "ERROR: multiple resume temp dirs found for ${out_4d}" >&2
      printf '%s\n' "${tmp_candidates[@]}" >&2
      exit 1
    elif (( ${#tmp_candidates[@]} == 1 )); then
      tmp_dir=${tmp_candidates[0]}
      log "Resuming warped volumes from $(basename "${tmp_dir}")"
    else
      tmp_dir=$(mktemp -d "${out_4d}.tmp.XXXXXX")
    fi

    local driver_prep_dir
    driver_prep_dir=$(prepare_driver_volumes "${in_4d}")
    local out_list=()
    local i idx in_vol out_vol
    for ((i=0; i<nvol; i++)); do
      printf -v idx "%06d" "${i}"
      in_vol="${driver_prep_dir}/vol_${idx}.nii.gz"
      out_vol="${tmp_dir}/out_${idx}.nii.gz"
      if [[ ! -f ${in_vol} ]]; then
        echo "ERROR: prepped volume missing: ${in_vol}" >&2
        exit 1
      fi
      if [[ ! -f ${out_vol} ]]; then
        antsApplyTransforms -d 3 -i "${in_vol}" -r "${ref_img}" -o "${out_vol}" -n Linear "${xforms[@]}"
      fi
      out_list+=("${out_vol}")
    done

    if (( ${#out_list[@]} == 0 )); then
      echo "ERROR: no warped fMRI volumes produced for ${in_4d}" >&2
      exit 1
    fi

    local scratch_parent=${SLURM_TMPDIR:-/tmp/${USER:-user}}
    mkdir -p "${scratch_parent}"
    local scratch_dir
    scratch_dir=$(mktemp -d "${scratch_parent%/}/concat_nifti_4d_${MOUSE_ID:-mouse}_XXXXXX")
    trap 'rm -rf "${scratch_dir}"' RETURN

    local local_output="${scratch_dir}/$(basename "${out_4d}")"
    local shared_copy_tmp="${out_4d}.scratchcopy.$$"

    # Rebuild the stage-local memmap from this output's warped volumes.
    # Do not reuse sibling concat_nifti_4d_*/stack.dat files here: those can
    # belong to a different reference space (for example struct vs Allen).
    python "${SCRIPT_DIR}/concat_nifti_4d.py" \
      --output "${local_output}" \
      --template-4d "${in_4d}" \
      --input-dir "${tmp_dir}" \
      --output-dtype float32

    cp "${local_output}" "${shared_copy_tmp}"
    mv "${shared_copy_tmp}" "${out_4d}"

    rm -rf "${tmp_dir}"
    rm -rf "${scratch_dir}"
    trap - RETURN
  }

  # source -> structural
  local struct_xforms=()
  if [[ -n ${mod_warp} ]]; then
    struct_xforms+=(-t "${mod_warp}")
  fi
  struct_xforms+=(-t "${mod_affine}")

  # source -> within-template
  local tpl_xforms=(-t "${struct_warp}" -t "${struct_affine}")
  if [[ -n ${mod_warp} ]]; then
    tpl_xforms+=(-t "${mod_warp}")
  fi
  tpl_xforms+=(-t "${mod_affine}")

  # source -> Allen
  local allen_xforms=(-t "${TPL2ALLEN_WARP}" -t "${TPL2ALLEN_AFFINE}" -t "${struct_warp}" -t "${struct_affine}")
  if [[ -n ${mod_warp} ]]; then
    allen_xforms+=(-t "${mod_warp}")
  fi
  allen_xforms+=(-t "${mod_affine}")

  if [[ ${src_ndim} -gt 1 ]]; then
    apply_4d_lowmem "${src}" "${struct_native}" "${out_in_struct}" "${struct_xforms[@]}"
    apply_4d_lowmem "${src}" "${WITHIN_TEMPLATE}" "${out_in_template}" "${tpl_xforms[@]}"
    apply_4d_lowmem "${src}" "${ATLAS_TEMPLATE}" "${out_in_allen}" "${allen_xforms[@]}"
  else
    apply_3d_once "${src}" "${struct_native}" "${out_in_struct}" "${struct_xforms[@]}"
    apply_3d_once "${src}" "${WITHIN_TEMPLATE}" "${out_in_template}" "${tpl_xforms[@]}"
    apply_3d_once "${src}" "${ATLAS_TEMPLATE}" "${out_in_allen}" "${allen_xforms[@]}"
  fi
}

sanitize_label() {
  local v=$1
  echo "${v}" | sed 's/[^A-Za-z0-9._-]/_/g'
}

run_one_image() {
  local tp=$1
  local role=$2
  local src=$3

  local base label out_dir
  base=$(basename "${src}")
  label=$(sanitize_label "${base%.nii.gz}")
  out_dir="${ALIGN_DIR}/${tp}/${role}_${label}"
  mkdir -p "${out_dir}"

  local ref3d="${out_dir}/${label}_ref3d.nii.gz"
  local pref="${out_dir}/${label}_to_struct_"
  local ref_in_struct="${out_dir}/${label}_ref3d_in_struct.nii.gz"
  local reg_log="${out_dir}/${label}_to_struct_registration.log"

  local src_ndim
  src_ndim=$(prepare_ref3d "${src}" "${ref3d}")
  register_modality_to_struct "${STRUCT_NATIVE_CUR}" "${ref3d}" "${pref}" "${ref_in_struct}" "${reg_log}"

  local mod_warp="${pref}1Warp.nii.gz"
  local mod_inv="${pref}1InverseWarp.nii.gz"
  local mod_aff="${pref}0GenericAffine.mat"
  [[ -f ${mod_aff} ]] || { echo "ERROR: missing ${mod_aff}" >&2; exit 1; }
  if [[ ! -f ${mod_warp} ]]; then
    mod_warp=""
    mod_inv=""
  fi

  local in_struct="${out_dir}/${label}_in_struct.nii.gz"
  local in_template="${out_dir}/${label}_in_template.nii.gz"
  local in_allen="${out_dir}/${label}_in_allen.nii.gz"
  if [[ ${src_ndim} -gt 1 ]]; then
    in_struct="${out_dir}/${label}_in_struct_4d.nii.gz"
    in_template="${out_dir}/${label}_in_template_4d.nii.gz"
    in_allen="${out_dir}/${label}_in_allen_4d.nii.gz"
  fi

  apply_chain "${src}" "${src_ndim}" "${STRUCT_NATIVE_CUR}" "${in_struct}" "${in_template}" "${in_allen}" \
    "${mod_warp}" "${mod_aff}" "${STRUCT_WARP_CUR}" "${STRUCT_AFF_CUR}"

  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "${tp}" "${role}" "${label}" "${src}" "${src_ndim}" "${ref3d}" "${STRUCT_NATIVE_CUR}" \
    "${in_struct}" "${in_template}" "${in_allen}" "${reg_log}" >>"${MANIFEST_TSV}"

  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "${tp}" "${role}" "${label}" "${mod_warp}" "${mod_inv}" "${mod_aff}" \
    "${STRUCT_WARP_CUR}" "${STRUCT_AFF_CUR}" "${TPL2ALLEN_WARP}" "${TPL2ALLEN_AFFINE}" >>"${CHAIN_TSV}"
}

mapfile -t SORTED_TPS < <(printf '%s\n' "${!STRUCT_NATIVE[@]}" | sort)
if [[ ${#SORTED_TPS[@]} -eq 0 ]]; then
  echo "ERROR: no timepoints parsed from anatomical transform manifest" >&2
  exit 1
fi

selected_images=0
for tp in "${SORTED_TPS[@]}"; do
  tp_dir="${MOUSE_ROOT}/${tp}"
  STRUCT_NATIVE_CUR="${STRUCT_NATIVE[${tp}]}"
  STRUCT_WARP_CUR="${STRUCT_TO_TEMPLATE_WARP[${tp}]}"
  STRUCT_AFF_CUR="${STRUCT_TO_TEMPLATE_AFFINE[${tp}]}"

  if [[ ! -f ${STRUCT_NATIVE_CUR} || ! -f ${STRUCT_WARP_CUR} || ! -f ${STRUCT_AFF_CUR} ]]; then
    log "WARN: ${tp}: missing structural chain file(s); skipping timepoint"
    continue
  fi
  if [[ ! -d ${tp_dir} ]]; then
    log "WARN: ${tp}: no matching timepoint directory at ${tp_dir}; skipping"
    continue
  fi

  t2=$(pick_t2_driver "${tp_dir}")
  mapfile -t fmri < <(
    {
      find "${tp_dir}" -maxdepth 1 -type f -name '*_fMRI_6000_*.nii.gz'
      find "${tp_dir}" -maxdepth 1 -type f -name '*_fMRI_7000_*.nii.gz'
    } | sort
  )
  fmri_joined=""
  if [[ ${#fmri[@]} -gt 0 ]]; then
    fmri_joined=$(printf '%s;' "${fmri[@]}")
    fmri_joined=${fmri_joined%;}
  fi
  printf '%s\t%s\t%s\t%s\t%s\n' "${tp}" "${STRUCT_NATIVE_CUR}" "${t2}" "${#fmri[@]}" "${fmri_joined}" >>"${SELECTED_TSV}"

  if [[ -n ${t2} ]]; then
    log "[${tp}] Aligning T2tse via anatomical chain: $(basename "${t2}")"
    run_one_image "${tp}" "t2tse" "${t2}"
    selected_images=$((selected_images + 1))
  else
    log "[${tp}] WARN: no T2tse scan found"
  fi

  if [[ ${#fmri[@]} -eq 0 ]]; then
    log "[${tp}] WARN: no fMRI 6000/7000 scans found"
  else
    for fmri_img in "${fmri[@]}"; do
      log "[${tp}] Aligning fMRI via anatomical chain: $(basename "${fmri_img}")"
      run_one_image "${tp}" "fmri" "${fmri_img}"
      selected_images=$((selected_images + 1))
    done
  fi
done

if [[ ${selected_images} -eq 0 ]]; then
  echo "ERROR: no T2tse/fMRI images were aligned (check availability and timepoint naming)" >&2
  exit 1
fi

log "Computing QC gates + fMRI run metrics"
python "${SCRIPT_DIR}/mouse_multimodal_qc.py" \
  --alignment-manifest "${MANIFEST_TSV}" \
  --out-dir "${TABLE_DIR}"

log "Done. Integrated multimodal outputs: ${OUT_MOUSE}"
log "Alignment manifest: ${MANIFEST_TSV}"
