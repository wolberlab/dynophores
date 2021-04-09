"""
dynophores.tests.fixures

Fixures to be used in unit testing.
"""

from pathlib import Path

import pytest

from dynophores import Dynophore

PATH_TEST_DATA = Path(__name__).parent / "dynophores/tests/data"


@pytest.fixture(scope="module")
def dynophore():

    dynophore = Dynophore.from_dir(PATH_TEST_DATA / "out")
    return dynophore


@pytest.fixture(scope="module")
def superfeature(dynophore):

    superfeature = dynophore.superfeatures["H[4615,4623,4622,4613,4621,4614]"]
    return superfeature


@pytest.fixture(scope="module")
def chemicalfeaturecloud3d(superfeature):

    cloud = superfeature.cloud
    return cloud


@pytest.fixture(scope="module")
def chemicalfeaturecloud3dpoint(superfeature):

    point = next(iter(superfeature.cloud.points))
    return point


@pytest.fixture(scope="module")
def envpartner(superfeature):

    envpartner = superfeature.envpartners["ILE-10-A[169,171,172]"]
    return envpartner
