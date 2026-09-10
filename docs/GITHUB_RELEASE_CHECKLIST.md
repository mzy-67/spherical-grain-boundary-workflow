# GitHub release checklist

Before the first public push:

1. Replace `OWNER` in `CITATION.cff` with the GitHub organization or username.
2. Confirm the repository name is `llzo-spherical-gb-workflow`.
3. Confirm contributor names and order with all coauthors.
4. Confirm the DeePMD checkpoint redistribution decision; it is excluded by default.
5. Confirm the Nernst-Einstein effective-volume expression against the final manuscript.
6. Add the remaining publication analysis scripts under `analysis/`.
7. Run `python scripts/validate_repository.py` in the `dp` environment.
8. Create the initial commit and push to GitHub.
9. Enable the GitHub repository in Zenodo and create release `v1.0.0` only after the reproducibility package is final.
10. Add the resulting software and dataset DOIs to `README.md`, `CITATION.cff`, and the manuscript.
