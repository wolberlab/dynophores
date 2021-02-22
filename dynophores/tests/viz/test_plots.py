"""
Unit tests for dynophore.viz.plots.

Will only test if plotting raises errors.
"""

import pytest
import matplotlib

from dynophores.viz import plot


@pytest.mark.parametrize(
    "superfeature_names, color_by_feature_type, max_frames",
    [
        ("all", True, 1000),  # Defaults
        ("all", False, 1000),
        (("all",), False, 1000),
        ("AR[4605,4607,4603,4606,4604]", False, 1000),
        (["AR[4605,4607,4603,4606,4604]", "AR[4622,4615,4623,4613,4614,4621]"], False, 1000),
        (["all", "AR[4622,4615,4623,4613,4614,4621]"], False, 1000),
    ],
)
def test_superfeatures_occurrences(
    dynophore, superfeature_names, color_by_feature_type, max_frames
):

    fig, ax = plot.superfeatures_occurrences(
        dynophore, superfeature_names, color_by_feature_type, max_frames
    )

    assert isinstance(fig, matplotlib.figure.Figure)
    assert isinstance(ax, matplotlib.axes.Subplot)

    func = plot.superfeatures_occurrences_interactive(dynophore)
    assert func is not None


@pytest.mark.parametrize(
    "superfeature_name, max_frames",
    [
        ("AR[4605,4607,4603,4606,4604]", 1000),
    ],
)
def test_envpartners_occurrences(dynophore, superfeature_name, max_frames):

    fig, ax = plot.envpartners_occurrences(dynophore, superfeature_name, max_frames)

    assert isinstance(fig, matplotlib.figure.Figure)
    assert isinstance(ax, matplotlib.axes.Subplot)


@pytest.mark.parametrize(
    "superfeature_name, max_frames",
    [
        ("AR[4605,4607,4603,4606,4604]", 1000),
    ],
)
def test_envpartners(dynophore, superfeature_name, max_frames):

    fig, axes = plot.envpartners(dynophore, superfeature_name, max_frames)
    print(axes.size)

    assert isinstance(fig, matplotlib.figure.Figure)
    assert axes.size == 4
    assert isinstance(axes[0][0], matplotlib.axes.Subplot)


def test_superfeatures_vs_envpartners(dynophore):

    fig, ax = plot.superfeatures_vs_envpartners(dynophore)

    assert isinstance(fig, matplotlib.figure.Figure)
    assert isinstance(ax, matplotlib.axes.Subplot)


@pytest.mark.parametrize(
    "superfeature_name, kind",
    [
        ("AR[4605,4607,4603,4606,4604]", "line"),
        ("AR[4605,4607,4603,4606,4604]", "hist"),
    ],
)
def test_envpartner_distances(dynophore, superfeature_name, kind):

    fig, ax = plot.envpartner_distances(dynophore, superfeature_name, kind)

    assert isinstance(fig, matplotlib.figure.Figure)
    assert isinstance(ax, matplotlib.axes.Subplot)


@pytest.mark.parametrize("superfeature_names", ["xxx"])
def test_superfeatures_occurrences_raises(dynophore, superfeature_names):

    with pytest.raises(KeyError):
        plot.superfeatures_occurrences(dynophore, superfeature_names)


@pytest.mark.parametrize("superfeature_name", ["xxx"])
def test_envpartners_occurrences_raises(dynophore, superfeature_name):

    with pytest.raises(KeyError):
        plot.envpartners_occurrences(dynophore, superfeature_name)


@pytest.mark.parametrize("superfeature_name", ["xxx"])
def test_envpartners_raises(dynophore, superfeature_name):

    with pytest.raises(KeyError):
        plot.envpartners(dynophore, superfeature_name)


@pytest.mark.parametrize(
    "superfeature_name, kind", [("xxx", "line"), ("AR[4605,4607,4603,4606,4604]", "xxx")]
)
def test_envpartner_distances_raises(dynophore, superfeature_name, kind):

    with pytest.raises(KeyError):
        plot.envpartner_distances(dynophore, superfeature_name, kind)
