# Mouse Longitudinal Workflow

Source guide: `Scripts/PYM_TECH/Parker/Wholeseome_Hero_Field_Guide.txt`

This bundle is a cleaned public copy of the Parker mouse workflow. It covers a
longitudinal pipeline designed to keep structural and functional mouse imaging
data on the same transformation chain across time.

At a high level, the project tries to solve a common imaging problem cleanly:
how to organize raw scan data, convert it reproducibly, build a within-subject
anatomical reference, and then reuse that reference to align other modalities
such as T2-weighted scans and resting-state fMRI.

## Scope

The workflow is organized around six stages:

1. Raw data extraction and inventory.
2. Subject/timepoint organization.
3. DICOM-to-NIfTI conversion.
4. Anatomical longitudinal alignment.
5. Integrated multimodal longitudinal processing.
6. Cohort-scale submission.

## How The Workflow Works

The pipeline starts with raw mouse imaging archives and turns them into a
subject-first longitudinal structure. DICOM metadata is parsed to understand
study chronology and scan content, then conversion scripts produce a stable
NIfTI layout. From there, one structural scan per timepoint is selected to
build a within-mouse anatomical template. That template becomes the common
reference frame for the rest of the workflow.

Once the anatomical chain exists, the multimodal scripts align T2-weighted
series and 4D fMRI runs into the same subject/template/atlas coordinate system.
QC scripts then summarize coverage, centroid behavior, and run-level metrics so
the output is inspectable at scale rather than only one mouse at a time.

## Representative Entry Points

These are good first files to read:

- `scripts/extract_nonhuman_mouse_dicom_metadata.py`
  metadata inventory and chronology extraction
- `scripts/organize_nonhuman_mouse_directory.py`
  subject/timepoint restructuring of raw scan folders
- `scripts/convert_nonhuman_mouse_organized_to_nifti.py`
  reproducible conversion into longitudinal NIfTI layout
- `scripts/select_mouse_driver_scans.py`
  selection of structural driver scans for the anatomical chain
- `scripts/mouse_longitudinal_pipeline.sh`
  core anatomical longitudinal workflow
- `scripts/mouse_longitudinal_fmri_t2tse.sh`
  multimodal integration using the established anatomical transform chain
- `scripts/mouse_multimodal_qc.py`
  QC gate and run-metric summarization

## What This Code Shows

- end-to-end workflow design from raw archive to aligned analysis-ready output
- explicit treatment of timepoint structure and multimodal registration
- mixing Python for structured data tasks with Bash for orchestration
- attention to operational concerns like batch execution, manifests, and
  low-memory handling for 4D data
- trying to make a complicated imaging workflow legible enough to rerun later

## Included Scripts

The `scripts/` folder contains 13 workspace-local files referenced by the
source guide, including:

- extraction and inventory scripts such as
  `run_unzip_mouse_mri_directory.sbatch`,
  `extract_nonhuman_mouse_dicom_metadata.py`, and
  `run_nonhuman_mouse_metadata.sbatch`
- organization and conversion tools such as
  `organize_nonhuman_mouse_directory.py`,
  `convert_nonhuman_mouse_organized_to_nifti.py`, and
  `run_convert_nonhuman_mouse_nifti.sbatch`
- anatomical workflow scripts such as `select_mouse_driver_scans.py`,
  `mouse_longitudinal_pipeline.sh`, `run_mouse_longitudinal_pipeline.sbatch`,
  and `mouse_longitudinal_metrics.py`
- multimodal processing and cohort submission scripts such as
  `mouse_longitudinal_fmri_t2tse.sh`, `mouse_multimodal_qc.py`, and
  `submit_mouse_longitudinal_cohort.sh`

## Not Bundled

These names appeared in the source guide but were not copied here:

- `QC_report.m`: not present as a standalone local script in this workspace
- `mouse_longitudinal_pipeline_README.m`: not present as a local script; the
  active workspace file is a markdown README rather than a MATLAB file

## Upload Notes

These copies were sanitized for public release. Workflow structure and tone
were preserved, while lab-specific paths and environment details were replaced
with placeholders where needed.
