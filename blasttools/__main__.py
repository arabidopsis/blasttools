from __future__ import annotations

from . import exact
from . import plants_ui
from .cli import blast

__all__ = [
    "plants_ui",
    "exact",
    "blast",
]

if __name__ == "__main__":
    blast.main(prog_name="blasttools")  # pylint: disable=no-value-for-parameter
