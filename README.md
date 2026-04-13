# Longitudinal Mouse Neuroimaging Pipeline

Longitudinal MRI in mice lets you track structural and functional brain changes
in the same animal over time — but only if the registration chain is consistent
enough that differences between timepoints reflect biology rather than
alignment error. This repository contains the pipeline I use to take raw mouse
MRI archives from scan day through anatomical template building, Allen atlas
registration, T2-weighted and 4D fMRI integration, and cohort-scale QC.

The pipeline is designed around one rule: every modality at every timepoint
shares the same anatomical transform chain. That keeps structural and
functional outputs in a single coordinate history so longitudinal cross-modality
comparisons do not drift apart.

<!-- SCREENSHOT 1 -->
<!-- Open a mouse brain image that has been registered to Allen atlas space.
     In Freeview or ITK-SNAP, load the registered structural image as the
     background and overlay the Allen atlas parcels (or just show the registered
     image on its own if parcels are not easily available).
     A coronal or axial slice showing the whole brain in atlas space works well.
     The goal is to show that the mouse brain is cleanly aligned to the atlas. -->
![Mouse brain in Allen atlas space](docs/img/mouse_allen_registration.png)
*Mouse structural image registered to Allen Mouse Brain Atlas space. All modalities and timepoints are ultimately brought into this common reference.*

---

## Why Longitudinal Mouse fMRI Is Tricky

Mice move between sessions. The scanner setup shifts. The brain itself changes
over the course of a study. Standard single-session registration approaches
fail here because they treat every timepoint independently, which means any
apparent longitudinal signal could just be registration inconsistency.

This pipeline builds a within-subject anatomical template from all timepoints
first, registers each timepoint to that template, and then registers the
template to the Allen atlas once. All modalities at all timepoints inherit that
chain. The 4D fMRI volumes are handled volume-by-volume through the same warp
to keep memory use manageable on large datasets.

<!-- SCREENSHOT 2 -->
<!-- In Freeview or ITK-SNAP, open structural images from two different
     timepoints for the same mouse, side by side or as two overlapping volumes
     with different color maps. The goal is to show that the anatomical
     alignment is consistent across time — the brain structures should line up.
     A coronal or axial slice at a recognizable landmark (hippocampus, striatum,
     cortex) works well.
     Screenshot showing both timepoints aligned. -->
![Longitudinal alignment across two timepoints](docs/img/longitudinal_alignment.png)
*Structural images from two timepoints for the same animal. Consistent alignment across sessions is the foundation for valid longitudinal comparison.*

---

## Requirements

- **ANTs** (2.4+) — anatomical template building and deformable registration
- **AFNI** — @SSwarper and supporting image operations
- **dcm2niix** — DICOM to NIfTI conversion
- **Python 3.8+** — `numpy`, `nibabel`, `scipy`, `pandas`
- Raw DICOM archives organized by subject and session
- Allen Mouse Brain Atlas template and parcel labels (50 µm or 25 µm)

**Atlas note:** The Allen atlas files are not included. Compatible files can be
downloaded from the Allen Brain Atlas data portal. The pipeline expects the
template at the path set by `ATLAS_TEMPLATE` and parcel labels at
`ATLAS_PARCELS`.

---

## Running the Pipeline

### Single mouse

```bash
# Anatomical longitudinal pipeline
Scripts/mouse_longitudinal_pipeline.sh \
  --mouse-root MRI_data/processed/nonhuman_longitudinal/<mouse_id> \
  --atlas-parcels <path_to_allen_parcels>

# Multimodal integration (runs after anatomical pipeline completes)
Scripts/mouse_longitudinal_fmri_t2tse.sh \
  --mouse-root MRI_data/processed/nonhuman_longitudinal/<mouse_id> \
  --anat-run-dir <path_to_completed_anatomical_output>
```

### Full cohort via SLURM

```bash
Scripts/submit_mouse_longitudinal_cohort.sh \
  --cohort-tsv MRI_data/processed/nonhuman_longitudinal/cohort_subjects.tsv \
  --atlas-parcels <path_to_allen_parcels>
```

---

## Pipeline Stages

### Stage 1 — Raw data extraction and inventory

`extract_nonhuman_mouse_dicom_metadata.py` reads DICOM headers and writes a
per-series and per-study metadata inventory. This runs before any conversion
and serves as a QC checkpoint on what was actually acquired.

### Stage 2 — Subject and timepoint organization

`organize_nonhuman_mouse_directory.py` restructures the raw DICOM tree into a
subject-first layout using symlinks and modality mapping by series prefix.
Outputs a manifest directly usable for downstream cohort submission.

### Stage 3 — DICOM to NIfTI conversion

`convert_nonhuman_mouse_organized_to_nifti.py` converts the organized tree to
a longitudinal NIfTI layout. Supports `--dry-run` and `--overwrite` and writes
a conversion manifest tracking the status of every series.

### Stage 4 — Anatomical longitudinal alignment

`mouse_longitudinal_pipeline.sh` is the core anatomical script. It selects one
structural driver scan per timepoint, builds a within-subject anatomical
template using ANTs, registers the template to the Allen atlas, and produces
Jacobian maps and parcel-level intensity tables for downstream analysis.

### Stage 5 — Multimodal integration

`mouse_longitudinal_fmri_t2tse.sh` integrates T2-weighted and 4D fMRI data
across timepoints by reusing the completed anatomical transform chain. Each
modality volume is registered to its same-timepoint structural image, then the
existing structural-to-template and template-to-Allen transforms are composed.
4D fMRI warps are applied volume-by-volume and restacked for low-memory safety.

### Stage 6 — QC

`mouse_multimodal_qc.py` computes explicit pass/fail QC gates and fMRI run
metrics — coverage, centroid distance, FD proxy, DVARS, tSNR, and censor
percentages — for every modality and timepoint in the cohort.

<!-- SCREENSHOT 3 -->
<!-- Open the qc_gates.tsv output from a completed pipeline run in a spreadsheet,
     or display it in the terminal with `column -t qc_gates.tsv | head -20`.
     The columns to show are: subject, timepoint, modality, pass, reason (if failed).
     Alternatively, open the QC_report.md in a rendered view.
     Screenshot showing the pass/fail gate rows for a few subjects. -->
![QC gate output](docs/img/qc_gates.png)
*QC gate output from `mouse_multimodal_qc.py`: per-modality pass/fail with coverage and centroid distance metrics.*

---

## Repository Layout

```
parker_mouse_longitudinal/    End-to-end longitudinal mouse pipeline
  scripts/                    Shell and Python scripts for all pipeline stages
  README.md                   Full workflow guide with stage-by-stage detail
```

---

## Scientific Context

This pipeline supports a longitudinal study of brain structure and function in
a mouse model of cardiac disease. The goal is to characterize how brain changes
track disease progression over time, using both structural volume changes
(Jacobians) and resting-state fMRI connectivity. Bringing T2-weighted and fMRI
data into the same anatomical coordinate history across timepoints is what makes
that comparison valid.

---

## Development History

This repository reflects the working pipeline as it was built and used. Commit
history tracks real decisions — including dead ends and corrections — not a
post-hoc reconstruction.

---

## Reproducibility Notes

Path conventions follow the internal study structure. Subject lists, atlas
files, and raw data are not redistributed. Adapting this pipeline to a new
dataset primarily means updating path and environment configuration; every
stage is documented in `parker_mouse_longitudinal/README.md`.
