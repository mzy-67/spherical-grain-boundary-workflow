# GitHub release checklist

Before the first public release:

1. Confirm the repository URL in `CITATION.cff` is `https://github.com/mzy-67/spherical-grain-boundary-workflow`.
2. Confirm the repository name and capitalization are final.
3. Confirm contributor names and order with all coauthors.
4. Confirm the DeePMD checkpoint redistribution decision; it is excluded by default.
5. Confirm the Nernst-Einstein region remains a radius-35 Å sphere intersected by a centered 30 Å-thick slab, as stated in the manuscript.
6. Run `python scripts/validate_repository.py` in the `dp` environment.
7. Create the initial commit and push to GitHub.
8. Make the repository public and create a tagged GitHub release.
9. Add the tagged software/data release URL to `README.md`, `CITATION.cff`, and the manuscript data-availability statement.
10. Remove any temporary write-enabled GitHub Deploy Key.
