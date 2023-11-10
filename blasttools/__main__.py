from __future__ import annotations

from . import exact  # noqa: pylint: disable=unused-import
from . import plants_ui  # noqa: pylint: disable=unused-import
from .cli import blast

if __name__ == "__main__":
    blast.main(prog_name="blasttools")  # pylint: disable=no-value-for-parameter
