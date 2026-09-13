# Changelog

## 0.2.0 - 2026-09-11

- Generalized chemistry, potential, carrier, charge, radii, and site-count settings.
- Renamed the primary package to `spherical_gb`; retained `jbfuncs` compatibility imports.
- Fixed antiparallel-vector rotation and enabled geometry/configuration tests in CI.
- Added scalar selected-radius outputs while preserving radial profiles.
- Retained LLZO as a documented reproduction example.

## Unreleased

- Added the 1,199-row property workbook and documented GitHub release assets containing 1,230 optimized/unoptimized structure pairs.
- Made the LAMMPS executable and manuscript MD protocol configurable.
- Added fail-fast handling for unsuccessful LAMMPS runs.
- Documented and tested the sphere–slab volume used for the Nernst-Einstein conversion.
- Added lightweight continuous-integration checks for source syntax, inputs, and publication metadata.
- Updated the GitHub data-availability documentation for `outputs.xlsx` and the released structure archives.

## 0.1.0 - 2026-09-10

- Created a publication-ready repository layout.
- Added the spherical GB Jobflow implementation and an unmodified provenance snapshot.
- Added the 1,291-entry GB orientation table and bulk reference inputs.
- Added configuration-driven HPC submission, environment metadata, citation metadata, and validation tests.
