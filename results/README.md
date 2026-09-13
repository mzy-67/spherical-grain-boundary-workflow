# Released results

This directory and the associated GitHub release contain the LLZO spherical grain-boundary results used by this project.

## Files

### `outputs.xlsx`

- One worksheet named `Sheet1`.
- 1 header row, 1,199 data rows, and 18 columns.
- GB IDs are unique and range from 1 to 1,291 with 92 sampled IDs absent because no complete property row was released.
- No empty cells occur in the 1,199 released rows.
- Column definitions and units are in [`../docs/DATA_DICTIONARY.md`](../docs/DATA_DICTIONARY.md).

### Structure archives

The large archives are assets of the [`dataset-v1.0` GitHub release](https://github.com/mzy-67/spherical-grain-boundary-workflow/releases/tag/dataset-v1.0):

- `opt-data1.zip` contains 1,230 optimized structures under `opt-data1/<GB id>.data`.
- `no-opt-data1.zip` contains the corresponding 1,230 unoptimized structures under `no-opt-data1/<GB id>.data`.
- The two archives contain identical GB-ID sets. They cover 1,230 of the 1,291 sampled orientations; 61 IDs have no released structure pair.
- All 1,199 IDs in `outputs.xlsx` occur in both archives. An additional 31 structure pairs have no complete row in `outputs.xlsx`.
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
b0079b5fabc0450e447bae502dc08ec5fab257c9683f18cd4f0ac3a33a9806bb  outputs.xlsx
a759fef19caa2f69bfd861cbc34d4b82ce02ecaace5b701cc9bf844568fcbfb2  opt-data1.zip
808bb010027f4adf326041cd3857db5351a1e27b002c6efe6b571714f7bd4bcf  no-opt-data1.zip
```

## Relationship between records

The numeric filename is the stable `id` used in `inputs/gb_orientations.csv` and `outputs.xlsx`. Join the files and table on this value. Missing IDs must not be interpreted as zero-valued properties or structures.

## Reuse notes

The conductivity columns are Nernst–Einstein estimates derived from the analyzed Li population. The 298.15 K value is an Arrhenius extrapolation from the simulated high-temperature values. Users should cite the associated article and the archived dataset/software DOI when available.
