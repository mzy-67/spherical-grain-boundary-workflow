"""Backward-compatible imports for the former :mod:`jbfuncs` package."""

from spherical_gb import SphericalGBWorkflowMaker

# Kept so existing LLZO scripts continue to run after the package rename.
SpheregbBOMaker = SphericalGBWorkflowMaker

__all__ = ["SphericalGBWorkflowMaker", "SpheregbBOMaker"]
__version__ = "0.2.0"
