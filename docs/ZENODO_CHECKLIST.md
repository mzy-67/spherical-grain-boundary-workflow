# Zenodo dataset checklist

Deposit the data as a separate Zenodo Dataset record and link its DOI from the GitHub software release. The current draft contains:

1. `outputs.xlsx`: one worksheet with 1,199 complete, unique GB records and 18 columns.
2. MD5 checksum `80905afe971f3753f9af15bf506c3cdf` (231,906 bytes), matching the local source file at upload time.

The full sampling table contains 1,291 candidates. The data record therefore represents the 1,199 candidates with complete released property rows; 92 sampled IDs are absent. See `docs/DATA_DICTIONARY.md` for field definitions and units.

Before publication:

1. Verify that every quantitative claim and figure in the article is reproducible from this record or from an explicitly linked companion record.
2. Include the software repository URL and, after the GitHub release is archived, its software DOI as related identifiers.
3. Check author names, affiliations, funding, keywords, access rights, and license in the Zenodo form.
4. Download the completed draft once and verify workbook integrity and checksums before clicking **Publish**.

Do not upload credentials, cluster configuration, redundant scheduler logs, or third-party potential files without redistribution permission.
