"""
Unit tests for dynophore.core.dynophore.Dynophore class.

Uses fixture tests.conftest.dynophore.
"""

from pathlib import Path

import pytest

from dynophores import Dynophore, SuperFeature

PATH_TEST_DATA = Path(__name__).parent / "dynophores" / "tests" / "data"


class TestsDynophore:
    """
    Test Dynophore class methods.
    """

    @pytest.mark.parametrize("id", ["1KE7-1"])
    def test_attributes(self, dynophore, id):

        assert dynophore.id == id
        assert isinstance(dynophore.superfeatures, list)
        for superfeature in dynophore.superfeatures:
            assert isinstance(superfeature, SuperFeature)

    @pytest.mark.parametrize("filepath", [PATH_TEST_DATA / "1KE7-1/DynophoreApp/data"])
    def test_from_files(self, filepath):

        dynophore = Dynophore.from_files(filepath)
        assert isinstance(dynophore, Dynophore)

    @pytest.mark.parametrize(
        "column_names, counts_sum",
        [
            (
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
                3141,
            )
        ],
    )
    def test_superfeatures_occurrences(self, dynophore, column_names, counts_sum):

        assert dynophore.superfeatures_occurrences.columns.to_list() == column_names
        assert dynophore.superfeatures_occurrences.index.to_list() == list(
            range(0, dynophore.n_frames)
        )
        assert dynophore.superfeatures_occurrences.sum().sum() == counts_sum

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

        assert list(dynophore.envpartners_occurrences.keys()) == list(counts_sum_dict.keys())
        counts_sum_dict_calculated = {
            superfeature: occurrences.sum().sum()
            for superfeature, occurrences in dynophore.envpartners_occurrences.items()
        }
        for (_, counts_sum_calculated), (_, counts_sum) in zip(
            counts_sum_dict_calculated.items(), counts_sum_dict.items()
        ):
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

        assert list(dynophore.envpartners_distances.keys()) == list(distances_sum_dict.keys())
        distances_sum_dict_calculated = {
            superfeature: occurrences.sum().sum()
            for superfeature, occurrences in dynophore.envpartners_distances.items()
        }
        for (_, distances_sum_calculated), (_, distances_sum) in zip(
            distances_sum_dict_calculated.items(), distances_sum_dict.items()
        ):
            assert pytest.approx(distances_sum_calculated) == distances_sum

    @pytest.mark.parametrize("data_type", ["xxx"])
    def test_envpartners_data_raises(self, dynophore, data_type):

        with pytest.raises(KeyError):
            dynophore._envpartners_data(data_type)

    @pytest.mark.parametrize("n_superfeatures", [10])
    def test_n_superfeatures(self, dynophore, n_superfeatures):

        assert dynophore.n_superfeatures == n_superfeatures

    @pytest.mark.parametrize("n_frames", [1002])
    def test_n_frames(self, dynophore, n_frames):

        assert dynophore.n_frames == n_frames

    @pytest.mark.parametrize(
        "count_sum, frequency_sum, superfeature_names, envpartner_names",
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
        self, dynophore, count_sum, frequency_sum, superfeature_names, envpartner_names
    ):

        # Count
        assert dynophore.count.columns.to_list() == superfeature_names
        assert dynophore.count.index.to_list() == envpartner_names
        assert dynophore.count.sum().sum() == count_sum

        # Frequency
        assert dynophore.frequency.columns.to_list() == superfeature_names
        assert dynophore.frequency.index.to_list() == envpartner_names
        assert pytest.approx(dynophore.frequency.sum().sum()) == frequency_sum

    @pytest.mark.parametrize(
        "valid_superfeature, superfeature_name",
        [(True, "AR[4605,4607,4603,4606,4604]"), (False, "xxx")],
    )
    def test_raise_keyerror_if_invalid_superfeature_name(
        self, dynophore, valid_superfeature, superfeature_name
    ):

        if valid_superfeature:
            assert dynophore.raise_keyerror_if_invalid_superfeature_name(superfeature_name) is None
        else:
            with pytest.raises(KeyError):
                dynophore.raise_keyerror_if_invalid_superfeature_name(superfeature_name)

    @pytest.mark.parametrize(
        "filepath, file_components",
        [
            (
                PATH_TEST_DATA
                / "1KE7-1/DynophoreApp/data"
                / "1KE7-1_data_superfeature_H[4599,4602,4601,4608,4609,4600]_100.0.txt",
                {
                    "filepath": PATH_TEST_DATA
                    / "1KE7-1/DynophoreApp/data"
                    / "1KE7-1_data_superfeature_H[4599,4602,4601,4608,4609,4600]_100.0.txt",
                    "dynophore_id": "1KE7-1",
                    "superfeature_id": "H[4599,4602,4601,4608,4609,4600]",
                    "superfeature_feature_type": "H",
                    "superfeature_atom_numbers": [4599, 4602, 4601, 4608, 4609, 4600],
                    "envpartner_id": None,
                    "envpartner_residue_name": None,
                    "envpartner_residue_number": None,
                    "envpartner_chain": None,
                    "envpartner_atom_numbers": None,
                },
            ),
            (
                PATH_TEST_DATA
                / "1KE7-1/DynophoreApp/data"
                / "1KE7-1_data_superfeature_HBA[4619]_12.3_envpartner_LYS_20_A[308]_4.1.txt",
                {
                    "filepath": PATH_TEST_DATA
                    / "1KE7-1/DynophoreApp/data"
                    / "1KE7-1_data_superfeature_HBA[4619]_12.3_envpartner_LYS_20_A[308]_4.1.txt",
                    "dynophore_id": "1KE7-1",
                    "superfeature_id": "HBA[4619]",
                    "superfeature_feature_type": "HBA",
                    "superfeature_atom_numbers": [4619],
                    "envpartner_id": "LYS-20-A[308]",
                    "envpartner_residue_name": "LYS",
                    "envpartner_residue_number": 20,
                    "envpartner_chain": "A",
                    "envpartner_atom_numbers": [308],
                },
            ),
        ],
    )
    def test_file_components(self, filepath, file_components):

        dynophore = Dynophore()
        file_components_calculated = dynophore._file_components(filepath)

        for (_, component_calculated), (_, component) in zip(
            file_components_calculated.items(), file_components.items()
        ):
            assert component_calculated == component
