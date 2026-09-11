# Scalable Spherical Grain-Boundary Workflow for LLZO

This repository contains the Jobflow-based workflow used to construct, optimize, anneal, and simulate an ensemble of spherical grain-boundary (GB) models of cubic Li7La3Zr2O12 (LLZO).

The workflow supports the study **An Ensemble of Spherical Grain-Boundary Models Enables Scalable Atomistic Investigation of LLZO Grain Boundaries** by Ziyi Man, Yaoshu Xie, Zhanlin Li, Lu Jiang, and Tingzheng Hou.

## Workflow

```text
Bulk LLZO structure
        |
Macroscopic GB geometry (misorientation and GB-plane normal)
        |
Construction of two spherical hemispherical grains
        |
Five-dimensional Bayesian optimization
        |
Local simulated annealing and structural relaxation
        |
MD at 973.15, 1073.15, 1173.15, and 1273.15 K
        |
GB energy, MSD, diffusion, conductivity, and anisotropy analysis
```

The Bayesian search uses five independent variables: a displacement of grain 1 along the GB normal, three translations of grain 2, and the interfacial gap. The selected solution is recorded as seven Cartesian/interface parameters `[x1, y1, z1, x2, y2, z2, gap]`.

## Repository contents

- `src/jbfuncs/gbmaker2.py`: spherical GB construction and Jobflow workflow.
- `scripts/submit_gb_workflow.py`: configuration-driven batch submission script.
- `inputs/gb_orientations.csv`: 1,291 sampled LLZO GB geometries.
- `inputs/sorted_optimized_bulk_POSCAR`: bulk LLZO input structure.
- `inputs/bulk_energy_reference.lammpstrj`: bulk atom-resolved reference data used by the workflow.
- `configs/workflow.example.yaml`: portable configuration template.
- `legacy/gbmaker2_original.py`: pre-publication source snapshot retained for provenance; one commented user-specific path was sanitized.

Large trajectory data, optimized structures, per-ion descriptors, trained models, and figure source data should be deposited separately as a Zenodo dataset. They are intentionally not tracked in Git.

## Software environment used for the reported calculations

- Python 3.12.12
- NumPy 2.1.3
- pandas 2.3.3
- pymatgen 2025.10.7
- jobflow 0.2.1
- jobflow-remote 0.1.8
- qtoolkit 0.1.6
- scikit-optimize 0.10.2
- interfacemaster 1.1.7
- matplotlib 3.10.8
- DeePMD-kit 3.0.2
- LAMMPS 29 Aug 2024 Update 1

`environment-lock.yml` records the exact Linux HPC environment used for provenance. It includes locally installed project packages and is not intended as a portable installer; use `environment.yml` for a clean installation.

## Installation

```bash
conda env create -f environment.yml
conda activate llzo-gb
python -m pip install -e .
```

`jobflow-remote` and `qtoolkit` may require the installation method used by your HPC facility if they are not available from your configured package channels.

## Interatomic potential

The DeePMD checkpoint is not redistributed in this repository. Set `potential.checkpoint_file` in a local configuration file to an authorized copy of the Li-La-Zr-O model used in the study. Cite the original potential publication when using it. Do not commit model files unless their redistribution license permits it.

## Configure a run

Copy the public template to a local configuration:

```bash
cp configs/workflow.example.yaml configs/workflow.local.yaml
```

Edit the checkpoint path and Jobflow-remote settings. `workflow.local.yaml` is excluded by `.gitignore` because it may contain cluster-specific paths.

The example configuration also records the MD protocol used in the manuscript: 4 ps equilibration followed by 50 ps production at 973.15, 1073.15, 1173.15, and 1273.15 K. The conductivity analysis uses Li ions in the intersection of a radius-35 Å sphere and a centered 30 Å-thick slab. Its geometric volume is evaluated analytically as

```text
V = 2π(R²h − h³/3),  R = 35 Å, h = 15 Å.
```

The LAMMPS executable defaults to `lmp` and can be changed through `simulation.lammps_executable` (for example, to `lmp_mpi`). A non-zero LAMMPS exit code stops the workflow immediately and reports the final output lines.

Validate without submitting:

```bash
python scripts/submit_gb_workflow.py \
  --config configs/workflow.local.yaml \
  --start-row 0 --stop-row 1 --dry-run
```

Submit rows 1250 through 1264 (zero-based indexing; `stop-row` is exclusive):

```bash
python scripts/submit_gb_workflow.py \
  --config configs/workflow.local.yaml \
  --start-row 1250 --stop-row 1265
```

## Reproducibility and provenance

The public module contains two non-scientific corrections relative to the preserved source snapshot:

1. A helper now refers to its `structure` argument rather than an undefined global variable.
2. A commented user-specific output path was replaced with a portable example path.

A pre-publication source snapshot is retained in `legacy/gbmaker2_original.py`; only a commented user-specific path was sanitized. No scientific constants or production calculation settings were changed during repository packaging.

The effective volume used for the Nernst-Einstein conversion is the sphere–slab intersection defined above, matching the manuscript method. Before tagging the archival release, confirm that the released source is the exact revision used for the reported results.

## Data availability

The orientation table in this repository contains 1,291 sampled GB candidates. The companion Zenodo draft contains `outputs.xlsx`, a table of 1,199 unique GBs for which all 18 released geometry, energetic, and transport fields are present. The remaining 92 sampled candidates do not have a complete row in the released property table. The workbook contains no structure archives; column definitions and units are documented in `docs/DATA_DICTIONARY.md`.

The article's final Data Availability statement should cite the published Zenodo DOI. Additional derived data used for figures or statistical claims should be included in that record or deposited as a clearly linked companion record.

## Citation

Citation metadata are provided in `CITATION.cff`. After GitHub-Zenodo archiving, add the software DOI to both `CITATION.cff` and this README.

## License and attribution

The code is distributed under the MIT License. The original JobflowFunctions copyright notice is retained. See `CONTRIBUTORS.md` for authorship and provenance.
