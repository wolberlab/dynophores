"""
Unit tests for dynophore.core.superfeature.SuperFeature class.

Uses fixture tests.conftest.superfeature.
"""

import pytest
import pandas as pd


class TestsSuperFeature:
    """
    Test SuperFeature class methods.
    """

    def test_envpartners_occurrences(self, superfeature):

        data = superfeature.envpartners_occurrences
        assert isinstance(data, pd.DataFrame)
        assert data.index.to_list() == list(range(0, superfeature.occurrences.size))
        assert data.columns.to_list() == [envpartner.id for envpartner in superfeature.envpartners]
        assert data.dtypes.unique() == "int64"

    def test_envpartners_distances(self, superfeature):

        data = superfeature.envpartners_distances
        assert isinstance(data, pd.DataFrame)
        assert data.index.to_list() == list(range(0, superfeature.occurrences.size))
        assert data.columns.to_list() == [envpartner.id for envpartner in superfeature.envpartners]
        assert data.dtypes.unique() == "float64"

    @pytest.mark.parametrize("data_type", ["xxx"])
    def test_data_raises(self, superfeature, data_type):

        with pytest.raises(KeyError):
            superfeature._data(data_type)

    @pytest.mark.parametrize("n_frames", [3])
    def test_n_frames(self, superfeature, n_frames):

        assert superfeature.n_frames == n_frames

    @pytest.mark.parametrize("count", [([2, 1, 1])])
    def test_count(self, superfeature, count):

        envpartner_ids = [envpartner.id for envpartner in superfeature.envpartners]
        count = pd.Series(count, index=["any"] + envpartner_ids)
        assert all(superfeature.count == count)

    @pytest.mark.parametrize("frequency", [([2 / 3.0, 1 / 3.0, 1 / 3.0])])
    def test_frequency(self, superfeature, frequency):

        envpartner_ids = [envpartner.id for envpartner in superfeature.envpartners]
        frequency = pd.Series(frequency, index=["any"] + envpartner_ids)
        frequency = round(frequency * 100, 2)
        assert all(superfeature.frequency == frequency)
