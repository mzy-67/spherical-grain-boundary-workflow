# Input files

- `gb_orientations.csv` contains 1,291 sampled GB geometries. Angles are in degrees. Rotation-axis and plane-normal components are Cartesian direction components inherited from the source table.
- `sorted_optimized_bulk_POSCAR` is the bulk LLZO structure used by the workflow.
- `bulk_energy_reference.lammpstrj` contains the atom-resolved bulk reference quantities used in GB excess calculations.

The DeePMD checkpoint is intentionally not included. Configure its authorized local path in `configs/workflow.local.yaml`.
