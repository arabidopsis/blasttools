from .cli import blast
from . import plants_ui  # pylint: disable=unused-import
from . import exact # pylint: disable=unused-import

if __name__ == "__main__":
    blast.main(prog_name="blasttools")  # pylint: disable=no-value-for-parameter
