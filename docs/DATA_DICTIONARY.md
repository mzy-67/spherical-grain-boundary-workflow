# Data dictionary

## `inputs/gb_orientations.csv`

| Column | Meaning |
| --- | --- |
| `id` | Stable GB identifier |
| `sigma` | Coincidence-site-lattice Σ value where applicable |
| `misorientation_angle_deg` | Grain misorientation angle in degrees |
| `axis_h`, `axis_k`, `axis_l` | Rotation-axis components |
| `plane_h`, `plane_k`, `plane_l` | GB-plane normal components |
| `gb_type` | Macroscopic classification: tilt, asymmetric tilt, twist, mixed, or unclassified |

## `results/outputs.xlsx`

`outputs.xlsx` contains one worksheet (`Sheet1`) with a header row and 1,199 data rows. Each `id` is unique. Starting from 1,291 sampled candidates, Bayesian optimization succeeded for 1,230 and failed for 61; MSD calculations then succeeded for 1,199 and failed for 31. The workbook contains the 1,199 candidates that completed MSD analysis.

| Column | Meaning | Unit |
| --- | --- | --- |
| `id` | Grain-boundary identifier | dimensionless |
| `sigma` | Coincidence-site-lattice Σ value | dimensionless |
| `Theta` | Misorientation angle | degree |
| `axis_h`, `axis_k`, `axis_l` | Rotation-axis indices | dimensionless |
| `BP_h`, `BP_k`, `BP_l` | Grain-boundary-plane indices | dimensionless |
| `gb_energy` | Grain-boundary energy | J m⁻² |
| `wb` | Work of separation | J m⁻² |
| `excess` | Excess volume per unit area | Å |
| `sigma_T_k_973.15` | Ionic conductivity at 973.15 K | S m⁻¹ |
| `sigma_T_k_1073.15` | Ionic conductivity at 1073.15 K | S m⁻¹ |
| `sigma_T_k_1173.15` | Ionic conductivity at 1173.15 K | S m⁻¹ |
| `sigma_T_k_1273.15` | Ionic conductivity at 1273.15 K | S m⁻¹ |
| `Ea_eV` | Arrhenius migration activation energy | eV |
| `sigma_T_k_25` | Ionic conductivity extrapolated to 298.15 K | S m⁻¹ |

The conductivity values use the Nernst–Einstein relation and the analyzed Li population and sphere–slab intersection volume documented in the repository README and workflow source. Despite the legacy column name, `sigma_T_k_25` denotes an Arrhenius extrapolation to 298.15 K (25 °C), not a direct simulation.

## GitHub release structure archives

The `dataset-v1.0` release also provides `no-opt-data1.zip` and `opt-data1.zip`. Each archive contains 1,230 LAMMPS data files named by the stable GB `id`; the former contains structures before Bayesian optimization and the latter contains the corresponding structures after Bayesian optimization. The Excel `id` maps directly to the numeric filename in both archives. Details and checksums are recorded in `results/README.md`.
