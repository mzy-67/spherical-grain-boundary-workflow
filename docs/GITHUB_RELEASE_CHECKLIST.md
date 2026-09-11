# GitHub release checklist

Before the first public release:

1. Confirm the repository URL in `CITATION.cff` is `https://github.com/mzy-67/LLZO-spherical-GB`.
2. Confirm the repository name and capitalization are final.
3. Confirm contributor names and order with all coauthors.
4. Confirm the DeePMD checkpoint redistribution decision; it is excluded by default.
5. Confirm the Nernst-Einstein region remains a radius-35 Å sphere intersected by a centered 30 Å-thick slab, as stated in the manuscript.
6. Run `python scripts/validate_repository.py` in the `dp` environment.
7. Create the initial commit and push to GitHub.
8. Make the repository public, enable it in Zenodo, and create a GitHub release whose tag matches the package and citation version (currently `v0.1.0`).
9. Add the resulting software and dataset DOIs to `README.md`, `CITATION.cff`, and the manuscript.
10. Revoke any temporary Zenodo token and remove any temporary write-enabled GitHub Deploy Key.
