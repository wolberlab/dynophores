"""
Unit tests for dynophore.viz.plot.interactive.

Will only test if interactive plotting raises errors.
"""

from dynophores.viz import plot


def test_superfeatures_vs_envpartners(dynophore):
    plot.interactive.superfeatures_vs_envpartners(dynophore)


def test_superfeatures_occurrences(dynophore):
    plot.interactive.superfeatures_occurrences(dynophore)


def test_envpartners_occurrences(dynophore):
    plot.interactive.envpartners_occurrences(dynophore)


def test_envpartners_distances(dynophore):
    plot.interactive.envpartners_distances(dynophore)


def test_envpartners_all_in_one(dynophore):
    plot.interactive.envpartners_all_in_one(dynophore)
