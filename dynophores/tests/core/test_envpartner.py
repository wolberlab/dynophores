"""
Unit tests for dynophore.core.envpartner.EnvPartner class.
"""

import pytest


class TestsEnvPartner:
    """
    Test EnvPartner class methods.
    """

    @pytest.mark.parametrize("n_frames", [3])
    def test_n_frames(self, envpartner, n_frames):

        assert envpartner.n_frames == n_frames

    @pytest.mark.parametrize("count", [1])
    def test_count(self, envpartner, count):

        assert envpartner.count == count

    @pytest.mark.parametrize("frequency", [1 / 3.0])
    def test_frequency(self, envpartner, frequency):

        frequency = round(frequency * 100, 2)
        assert envpartner.frequency == frequency
