# Spherical Grain-Boundary Workflow

A material-agnostic, configuration-driven workflow for constructing, optimizing, annealing, and analyzing spherical grain-boundary (GB) models.

The code uses pymatgen and InterfaceMaster for geometry generation, scikit-optimize for interfacial translation/gap searches, Jobflow for orchestration, and LAMMPS for energy minimization and molecular dynamics. The interatomic potential, chemical species, masses, mobile carrier, charge number, temperature schedule, radii, and cluster resources are all configured in YAML rather than hard-coded in Python.

## Workflow

```text
Bulk crystal + orientation table
        |
Two spherical hemispherical grains
        |
Bayesian optimization of translation and gap
        |
Local annealing and structural relaxation
        |
Radial GB-energy / separation-work / excess-volume profiles
        |
Optional carrier MSD, diffusion, conductivity, and Arrhenius analysis
```

## What is generic now

- Any composition supported by pymatgen can be used through `material.species_order`.
- LAMMPS atom types, masses, and `dump_modify ... element` names are generated dynamically.
- The mobile species and its signed charge number are configurable.
- Potential commands are configurable through `potential.pair_style`, `pair_style_args`, and `pair_coeff`.
- Bulk-site indexing uses the actual number of sites in the input structure, not a fixed LLZO count.
- Energy radii are generated from the configured endpoint and include `grain_boundary.energy_radius`.
- The antiparallel-normal rotation case is handled correctly and tested.
- Annealing returns both radial profiles and explicit scalar values at the selected radius.

## Repository layout

- `src/spherical_gb/`: generic Python package.
- `src/jbfuncs/`: compatibility imports for scripts written against the original release.
- `scripts/submit_gb_workflow.py`: YAML-driven validation and Jobflow submission.
- `configs/workflow.example.yaml`: material-neutral template.
- `examples/llzo/`: settings that reproduce the original LLZO workflow.
- `examples/single_gb/`: one-orientation LLZO smoke test.
- `inputs/`: LLZO example/reference inputs retained for reproducibility.
- `legacy/`: original pre-publication source snapshot; not used by the package.

## Installation

```bash
conda env create -f environment.yml
conda activate spherical-gb
python -m pip install -e .
```

`jobflow-remote`, `qtoolkit`, LAMMPS, and any LAMMPS potential plugin (for example DeePMD-kit) must also be available on the target machine or HPC cluster.

## Configure a material

Copy the generic template:

```bash
cp configs/workflow.example.yaml configs/workflow.local.yaml
```

At minimum, set:

```yaml
inputs:
  orientation_file: /path/to/orientations.csv
  structure_file: /path/to/bulk_structure.cif
  bulk_energy_trajectory: /path/to/bulk_reference.lammpstrj

material:
  species_order: [Na, Cl]
  masses: {}                 # optional overrides; pymatgen supplies defaults
  mobile_species: Na
  charge_number: 1
  origin_anchor_species: null

potential:
  pair_style: deepmd
  checkpoint_file: /path/to/model.pb
  pair_style_args: "{checkpoint_file}"
  pair_coeff: "* *"
```

`species_order` defines LAMMPS types `1..N` and must exactly match the elements in the structure. For a non-DeePMD potential, set `pair_style`, `pair_style_args`, and `pair_coeff` to the corresponding LAMMPS command arguments. The `{checkpoint_file}` placeholder is optional.

The orientation CSV must contain:

```text
id,sigma,misorientation_angle_deg,axis_h,axis_k,axis_l,plane_h,plane_k,plane_l,gb_type
```

Validate one orientation without submission:

```bash
python scripts/submit_gb_workflow.py \
  --config configs/workflow.local.yaml \
  --start-row 0 --stop-row 1 --dry-run
```

Remove `--dry-run` to submit. Row indices are zero-based and `--stop-row` is exclusive.

## LLZO reproduction example

The original Li7La3Zr2O12 inputs and its 1,291-orientation table remain available as an example. Copy `examples/llzo/config.yaml`, set the authorized DeePMD checkpoint and your scheduler settings, then run the same submit command. The companion research dataset contains `outputs.xlsx` with 1,199 completed records; 92 runs failed during relaxation or MD and are not represented in that table.

The DeePMD checkpoint and large production outputs are not redistributed here. Cite the original potential publication and follow its license.

## Output contract

The annealing job returns:

- `radii_A` and three radial profiles;
- `selected_radius_A`;
- scalar `gb_energy`, `work_of_separation`, and `excess_volume` at that radius;
- legacy `wb` and `excess` arrays for older notebooks.

Transport output records the configured carrier count as `n_mobile`, its number density, charge number, diffusion coefficients, conductivity values, and Arrhenius fit.

## Limitations

A local test suite can verify configuration, geometry helpers, metadata, and source invariants. A complete scientific validation still requires a licensed potential, a compatible LAMMPS build, and representative HPC runs for the chosen material.

## Citation and license

Citation metadata are in `CITATION.cff`. The software is MIT licensed. The LLZO example supports the associated study *An Ensemble of Spherical Grain-Boundary Models Enables Scalable Atomistic Investigation of LLZO Grain Boundaries*.
