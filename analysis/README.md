# Analysis code

Place publication analysis in the following directories before the archival release:

- `gb_energy/`: GB-energy aggregation and structural descriptors.
- `transport/`: MSD, diffusion, conductivity, activation energy, and anisotropy.
- `energy_ml/`: machine-learning prediction of GB energy.
- `softness/`: mobility labeling, descriptor generation, linear-SVM training, and evaluation.
- `cavity_network/`: cavity construction, hopping-network mapping, and topology metrics.
- `figures/`: scripts that reproduce main and supporting figures from the Zenodo dataset.

Each script should accept explicit input/output arguments and avoid HPC-specific absolute paths.
