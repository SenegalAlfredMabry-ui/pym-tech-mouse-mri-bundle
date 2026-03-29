# Mouse MRI/fMRI Project Bundle

This repository is a sanitized code bundle for my longitudinal mouse MRI and
fMRI workflow. The code covers the pipeline from raw imaging archives through
organized subject/timepoint structure, DICOM-to-NIfTI conversion, anatomical
alignment, multimodal integration, QC, and cohort-scale submission.

The project is organized as a workflow bundle rather than a general-purpose
software package. It is meant to show how I structure imaging pipelines that
need to stay reproducible across multiple scan types and timepoints.

## Repository Layout

- `parker_mouse_longitudinal/`
  end-to-end longitudinal mouse pipeline, including raw-data inventory,
  organization, conversion, anatomical processing, multimodal alignment, QC,
  and batch submission

## What A Reviewer Can Look For

- pipeline decomposition into clean stages instead of one monolithic script
- mixed Bash and Python workflow design for real imaging operations
- deliberate handling of multimodal transform chains rather than isolated
  one-off registrations
- QC generation as a first-class output rather than an afterthought
- awareness of compute constraints, especially for 4D fMRI handling

## Public Release Notes

This is a public-safe copy of the original internal bundle.

- workflow structure and naming style were preserved
- lab-specific paths were replaced with placeholders
- environment-specific details were generalized where needed

This makes the repository suitable as a code sample and project overview, even
though some local configuration would still be required to run it directly.
