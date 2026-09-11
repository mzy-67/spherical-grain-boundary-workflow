"""Material-agnostic spherical grain-boundary workflow."""

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .workflow import SphericalGBWorkflowMaker

__all__ = ["SphericalGBWorkflowMaker"]
__version__ = "0.2.0"


def __getattr__(name: str):
    if name == "SphericalGBWorkflowMaker":
        from .workflow import SphericalGBWorkflowMaker

        return SphericalGBWorkflowMaker
    raise AttributeError(name)
