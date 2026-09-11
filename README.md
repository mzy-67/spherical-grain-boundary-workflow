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

Copy the example configuration:

```bash
cp configs/workflow.example.yaml configs/workflow.local.yaml
```

Set the DeePMD checkpoint and Jobflow-remote options in `workflow.local.yaml`. This local file is ignored by Git because it may contain cluster-specific paths. The LAMMPS executable can be changed with `simulation.lammps_executable`.

Validate without submitting:

```bash
python scripts/submit_gb_workflow.py \
  --config configs/workflow.local.yaml \
  --start-row 0 --stop-row 1 --dry-run
```

Remove `--dry-run` to submit the selected rows. Row indices are zero-based and `--stop-row` is exclusive.

The example configuration contains the manuscript MD settings: 4 ps equilibration and 50 ps production at 973.15, 1073.15, 1173.15, and 1273.15 K. Conductivity is evaluated for Li ions within the intersection of a radius-35 Å sphere and a centered 30 Å slab, with volume

```text
V = 2π(R²h − h³/3),  R = 35 Å, h = 15 Å.
```

## Data availability

The orientation table contains 1,291 sampled GB candidates. The companion Zenodo dataset contains `outputs.xlsx` with 1,199 completed geometry, energetic, and transport records. The other 92 calculations failed during structural relaxation or molecular dynamics and are not included. Column definitions and units are provided in [`docs/DATA_DICTIONARY.md`](docs/DATA_DICTIONARY.md).

The reserved dataset DOI is `10.5281/zenodo.22687083`; it will become active when the Zenodo record is published.

## Citation

Citation metadata are provided in [`CITATION.cff`](CITATION.cff).

## License and attribution

The code is distributed under the MIT License. The original JobflowFunctions copyright notice is retained. See `CONTRIBUTORS.md` for authorship and provenance.
