"""
Unit tests for dynophore.viz.view3d.parser.
"""

from pathlib import Path

import pytest

from dynophores.viz import view3d

PATH_TEST_DATA = Path(__name__).parent / "dynophores" / "tests" / "data"


@pytest.mark.parametrize(
    "filepath",
    [PATH_TEST_DATA / "1KE7-1/DynophoreApp/dynophore.pml"],
)
def test_parse_pml(filepath):
    dynophore_dict = view3d.parser._parse_pml(filepath)
    for superfeature_name, coordinates in dynophore_dict.items():
        # TODO will need more specific unit testing
        assert isinstance(superfeature_name, str)
        assert isinstance(coordinates, list)
        assert isinstance(coordinates[0], list)
        assert isinstance(coordinates[0][0], float)
