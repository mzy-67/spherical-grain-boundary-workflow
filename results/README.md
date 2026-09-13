# Released results

This directory and the associated GitHub release contain the LLZO spherical grain-boundary results used by this project.

## Files

### `outputs.xlsx`

- One worksheet named `Sheet1`.
- 1 header row, 1,211 data rows, and 16 columns.
- GB IDs are unique and range from 1 to 1,291. The workbook contains the 1,211 candidates that completed both Bayesian optimization and MSD calculation.
- No empty cells occur in the 1,211 released rows.
- Column definitions and units are in [`../docs/DATA_DICTIONARY.md`](../docs/DATA_DICTIONARY.md).

### Structure archives

The large archives are assets of the [`dataset-v1.0` GitHub release](https://github.com/mzy-67/spherical-grain-boundary-workflow/releases/tag/dataset-v1.0):

- `no-opt-data1.zip` contains 1,230 pre-Bayesian-optimization structures under `no-opt-data1/<GB id>.data`.
- `opt-data1.zip` contains the corresponding 1,230 post-Bayesian-optimization structures under `opt-data1/<GB id>.data`.
- The two archives contain identical GB-ID sets. Bayesian optimization succeeded for these 1,230 candidates and failed for 61 of the original 1,291 candidates.
- MSD calculations succeeded for 1,211 of the 1,230 optimized structures and failed for 19. These 1,211 successful MSD results form the rows of `outputs.xlsx`.
- Files use the LAMMPS data format with `units metal`. Atom types are ordered as Li, La, Zr, and O (`1` through `4`).

The structure archives are large (approximately 1.6 GB and 970 MB compressed) and are intentionally excluded from Git history.

## Download and verify

```bash
curl -LO https://github.com/mzy-67/spherical-grain-boundary-workflow/releases/download/dataset-v1.0/opt-data1.zip
curl -LO https://github.com/mzy-67/spherical-grain-boundary-workflow/releases/download/dataset-v1.0/no-opt-data1.zip
shasum -a 256 results/outputs.xlsx opt-data1.zip no-opt-data1.zip
```

Expected SHA-256 checksums:

```text
9bb1a4c87d5e662409e0b9f9533a776979e30de2699e7813846d1084e186fae0  outputs.xlsx
a759fef19caa2f69bfd861cbc34d4b82ce02ecaace5b701cc9bf844568fcbfb2  opt-data1.zip
808bb010027f4adf326041cd3857db5351a1e27b002c6efe6b571714f7bd4bcf  no-opt-data1.zip
```

## Relationship between records

The numeric filename is the stable `id` used in `inputs/gb_orientations.csv` and the `id` column of `outputs.xlsx`. For example, the Excel row with `id = 100` corresponds to `no-opt-data1/100.data` and `opt-data1/100.data`. Join the pre-optimization structure, post-optimization structure, and calculated result on this value. Missing IDs represent failed workflow stages and must not be interpreted as zero-valued properties or structures.

## Reuse notes

The conductivity columns are Nernst–Einstein estimates derived from the analyzed Li population. Users should cite the associated article and the tagged GitHub data/software release.
