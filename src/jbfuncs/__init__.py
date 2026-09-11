"""Jobflow utilities for spherical LLZO grain-boundary calculations."""

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .gbmaker2 import SpheregbBOMaker

__all__ = ["SpheregbBOMaker"]
__version__ = "0.1.0"


def __getattr__(name: str):
    if name == "SpheregbBOMaker":
        from .gbmaker2 import SpheregbBOMaker

        return SpheregbBOMaker
    raise AttributeError(name)
