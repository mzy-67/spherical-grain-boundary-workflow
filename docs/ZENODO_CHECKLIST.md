# Zenodo dataset checklist

Deposit the data as a separate Zenodo Dataset record and link its DOI from the GitHub software release. The GitHub repository/release now provides:

1. `results/outputs.xlsx`: one worksheet with 1,199 complete, unique GB records and 18 columns.
2. `opt-data1.zip`: 1,230 optimized LAMMPS structures.
3. `no-opt-data1.zip`: 1,230 corresponding unoptimized LAMMPS structures.

The authoritative SHA-256 checksums are recorded in `results/README.md`. Mirror all three files to the Zenodo dataset, or make the GitHub release an explicit related identifier if the structure archives are intentionally hosted only on GitHub.

The full sampling table contains 1,291 candidates. The data record therefore represents the 1,199 candidates with complete released property rows; 92 sampled IDs are absent. See `docs/DATA_DICTIONARY.md` for field definitions and units.

Before publication:

1. Verify that every quantitative claim and figure in the article is reproducible from this record or from an explicitly linked companion record.
2. Include the software repository URL and, after the GitHub release is archived, its software DOI as related identifiers.
3. Check author names, affiliations, funding, keywords, access rights, and license in the Zenodo form.
4. Download the completed draft once and verify workbook/archive integrity and checksums before clicking **Publish**.

Do not upload credentials, cluster configuration, redundant scheduler logs, or third-party potential files without redistribution permission.
