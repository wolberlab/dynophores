"""
Unit tests for dynophore.core.dynophore.Dynophore class.

Uses fixture tests.conftest.dynophore.
"""

from pathlib import Path

import pytest

from dynophores import parsers
from dynophores import Dynophore
from dynophores.core.superfeature import SuperFeature
from dynophores.core.ligand import Ligand


PATH_TEST_DATA = Path(__name__).parent / "dynophores/tests/data"


class TestsDynophore:
    """
    Test Dynophore class methods.
    """

    @pytest.mark.parametrize("id", ["dynophore_1KE7"])
    def test_attributes(self, dynophore, id):

        dynophore_dict = parsers._json_pml_to_dict(
            PATH_TEST_DATA / "out/1KE7_dynophore.json",
            PATH_TEST_DATA / "out/1KE7_dynophore.pml",
        )
        dynophore = Dynophore(**dynophore_dict)
        assert isinstance(dynophore, Dynophore)
        assert list(dynophore.__dict__) == ["id", "ligand", "superfeatures"]

        # Test class attributes - check for data types
        assert isinstance(dynophore.id, str)
        assert isinstance(dynophore.ligand, Ligand)
        assert list(dynophore.ligand.__dict__) == [
            "name",
            "smiles",
            "mdl_mol_buffer",
            "atom_serials",
        ]
        assert isinstance(dynophore.superfeatures, dict)
        assert isinstance(next(iter(dynophore.superfeatures.values())), SuperFeature)

    @pytest.mark.parametrize("filepath", [PATH_TEST_DATA / "out"])
    def test_from_dir(self, filepath):

        dynophore = Dynophore.from_dir(filepath)
        assert isinstance(dynophore, Dynophore)

    @pytest.mark.parametrize(
        "column_names, counts_sum",
        [
            (
                [
                    "H[4599,4602,4601,4608,4609,4600]",
                    "H[4615,4623,4622,4613,4621,4614]",
                    "HBA[4596]",
                    "HBA[4619]",
                    "HBD[4612]",
                    "AR[4622,4615,4623,4613,4614,4621]",
                    "HBA[4606]",
                    "HBD[4598]",
                    "HBA[4618]",
                    "AR[4605,4607,4603,4606,4604]",
                ],
                3141,
            )
        ],
    )
    def test_superfeatures_occurrences(self, dynophore, column_names, counts_sum):

        data = dynophore.superfeatures_occurrences
        assert data.columns.to_list() == column_names
        assert data.index.to_list() == list(range(0, dynophore.n_frames))
        assert data.sum().sum() == counts_sum
        assert data.dtypes.unique() == "int32"

    @pytest.mark.parametrize(
        "counts_sum_dict",
        [
            {
                "AR[4605,4607,4603,4606,4604]": 1,
                "AR[4622,4615,4623,4613,4614,4621]": 40,
                "HBA[4596]": 862,
                "HBA[4606]": 20,
                "HBA[4618]": 2,
                "HBA[4619]": 126,
                "HBD[4598]": 10,
                "HBD[4612]": 80,
                "H[4599,4602,4601,4608,4609,4600]": 5093,
                "H[4615,4623,4622,4613,4621,4614]": 1997,
            }
        ],
    )
    def test_envpartners_occurrences(self, dynophore, counts_sum_dict):

        assert sorted(list(dynophore.envpartners_occurrences.keys())) == sorted(
            list(counts_sum_dict.keys())
        )
        for superfeature, occurrences in dynophore.envpartners_occurrences.items():
            counts_sum_calculated = occurrences.sum().sum()
            counts_sum = counts_sum_dict[superfeature]
            assert counts_sum_calculated == counts_sum

    @pytest.mark.parametrize(
        "distances_sum_dict",
        [
            {
                "AR[4605,4607,4603,4606,4604]": 5405.736763,
                "AR[4622,4615,4623,4613,4614,4621]": 21501.949293,
                "HBA[4596]": 3352.240247,
                "HBA[4606]": 11580.393482,
                "HBA[4618]": 5917.65908,
                "HBA[4619]": 28431.781275,
                "HBD[4598]": 9999.710108,
                "HBD[4612]": 30733.293196,
                "H[4599,4602,4601,4608,4609,4600]": 30959.400883,
                "H[4615,4623,4622,4613,4621,4614]": 15190.159939,
            }
        ],
    )
    def test_envpartners_distances(self, dynophore, distances_sum_dict):

        assert sorted(list(dynophore.envpartners_distances.keys())) == sorted(
            list(distances_sum_dict.keys())
        )
        for superfeature, distances in dynophore.envpartners_distances.items():
            distances_sum_calculated = distances.sum().sum()
            distances_sum = distances_sum_dict[superfeature]
            assert pytest.approx(distances_sum_calculated) == distances_sum

    @pytest.mark.parametrize("n_superfeatures", [10])
    def test_n_superfeatures(self, dynophore, n_superfeatures):

        assert dynophore.n_superfeatures == n_superfeatures

    @pytest.mark.parametrize("n_frames", [1002])
    def test_n_frames(self, dynophore, n_frames):

        assert dynophore.n_frames == n_frames

    @pytest.mark.parametrize(
        "count_sum, frequency_sum, superfeature_ids, envpartner_names",
        [
            (
                11372,
                1134.96,
                [
                    "AR[4605,4607,4603,4606,4604]",
                    "AR[4622,4615,4623,4613,4614,4621]",
                    "HBA[4596]",
                    "HBA[4606]",
                    "HBA[4618]",
                    "HBA[4619]",
                    "HBD[4598]",
                    "HBD[4612]",
                    "H[4599,4602,4601,4608,4609,4600]",
                    "H[4615,4623,4622,4613,4621,4614]",
                ],
                [
                    "ALA-144-A[2263,2266]",
                    "ALA-31-A[488,491]",
                    "ASP-86-A[1313]",
                    "ASP-86-A[1319]",
                    "ASP-86-A[1320]",
                    "GLN-131-A[2057]",
                    "GLN-131-A[2061]",
                    "GLN-131-A[2062]",
                    "GLU-81-A[1228]",
                    "HIS-84-A[1284,1285,1286,1287,1288]",
                    "ILE-10-A[165]",
                    "ILE-10-A[169,171,172]",
                    "ILE-10-A[169,171]",
                    "LEU-134-A[2109,2110,2111]",
                    "LEU-83-A[1260]",
                    "LEU-83-A[1263]",
                    "LYS-129-A[2026]",
                    "LYS-20-A[308]",
                    "LYS-20-A[316]",
                    "LYS-89-A[1374]",
                    "PHE-82-A[1245,1246,1247,1248,1249,1250]",
                    "VAL-18-A[275,276,277]",
                    "any",
                ],
            )
        ],
    )
    def test_count_frequency(
        self, dynophore, count_sum, frequency_sum, superfeature_ids, envpartner_names
    ):

        # TODO remove this when fixed in DynophoreApp json export
        envpartner_names = [i.replace("_", "-") for i in envpartner_names]

        # Count
        assert sorted(dynophore.count.columns.to_list()) == superfeature_ids
        assert dynophore.count.index.to_list() == envpartner_names
        assert dynophore.count.sum().sum() == count_sum
        assert dynophore.count.dtypes.unique() == "int32"

        # Frequency
        assert sorted(dynophore.frequency.columns.to_list()) == superfeature_ids
        assert dynophore.frequency.index.to_list() == envpartner_names
        assert pytest.approx(dynophore.frequency.sum().sum()) == frequency_sum

    @pytest.mark.parametrize(
        "valid_superfeature, superfeature_id",
        [(True, "AR[4605,4607,4603,4606,4604]"), (False, "xxx")],
    )
    def test_raise_keyerror_if_invalid_superfeature_id(
        self, dynophore, valid_superfeature, superfeature_id
    ):

        if valid_superfeature:
            assert dynophore._raise_keyerror_if_invalid_superfeature_id(superfeature_id) is None
        else:
            with pytest.raises(KeyError):
                dynophore._raise_keyerror_if_invalid_superfeature_id(superfeature_id)

    def test_superfeatures_atom_serials(self, dynophore):

        atom_serials = dynophore.superfeatures_atom_serials
        assert isinstance(atom_serials, dict)
        assert list(atom_serials) == list(dynophore.superfeatures)
        assert isinstance(next(iter(atom_serials.values()))[0], int)

    def test_superfeatures_colors(self, dynophore):

        colors = dynophore.superfeatures_colors
        assert isinstance(colors, dict)
        assert list(colors) == list(dynophore.superfeatures)
        assert isinstance(next(iter(colors.values())), str)
        assert len(next(iter(colors.values()))) in [6, 7]  # With/without #-prefix
