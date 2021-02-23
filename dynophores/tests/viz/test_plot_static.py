"""
Unit tests for dynophore.viz.plot.static.

Will only test if static plotting raises errors.
"""

import pytest
import matplotlib

from dynophores.viz import plot


@pytest.mark.parametrize(
    "superfeature_names",
    [
        ("all"),  # Default
        (("all",)),
        ("AR[4605,4607,4603,4606,4604]"),
        (["AR[4605,4607,4603,4606,4604]", "AR[4622,4615,4623,4613,4614,4621]"]),
        (["all", "AR[4622,4615,4623,4613,4614,4621]"]),
    ],
)
def test_superfeatures_vs_envpartners(dynophore, superfeature_names):

    fig, ax = plot.static.superfeatures_vs_envpartners(dynophore, superfeature_names)

    assert isinstance(fig, matplotlib.figure.Figure)
    assert isinstance(ax, matplotlib.axes.Subplot)


@pytest.mark.parametrize("superfeature_names", ["xxx"])
def test_superfeatures_vs_envpartners_raises(dynophore, superfeature_names):

    with pytest.raises(KeyError):
        plot.static.superfeatures_vs_envpartners(dynophore, superfeature_names)


@pytest.mark.parametrize(
    "superfeature_names, color_by_feature_type, n_equidistant_frames",
    [
        ("all", True, 1000),  # Defaults
        ("all", False, 1000),
        (("all",), False, 1000),
        ("AR[4605,4607,4603,4606,4604]", False, 1000),
        # ("AR[4605,4607,4603,4606,4604]", False, 1),
        (["AR[4605,4607,4603,4606,4604]", "AR[4622,4615,4623,4613,4614,4621]"], False, 1000),
        (["all", "AR[4622,4615,4623,4613,4614,4621]"], False, 1000),
    ],
)
def test_superfeatures_occurrences(
    dynophore, superfeature_names, color_by_feature_type, n_equidistant_frames
):

    fig, ax = plot.static.superfeatures_occurrences(
        dynophore, superfeature_names, color_by_feature_type, n_equidistant_frames
    )

    assert isinstance(fig, matplotlib.figure.Figure)
    assert isinstance(ax, matplotlib.axes.Subplot)


@pytest.mark.parametrize("superfeature_names", ["xxx"])
def test_superfeatures_occurrences_raises(dynophore, superfeature_names):

    with pytest.raises(KeyError):
        plot.static.superfeatures_occurrences(dynophore, superfeature_names)


@pytest.mark.parametrize(
    "superfeature_names, n_equidistant_frames",
    [
        ("AR[4605,4607,4603,4606,4604]", 1000),
        (["AR[4605,4607,4603,4606,4604]", "AR[4622,4615,4623,4613,4614,4621]"], 1000),
    ],
)
def test_envpartners_occurrences(dynophore, superfeature_names, n_equidistant_frames):

    fig, axes = plot.static.envpartners_occurrences(
        dynophore, superfeature_names, n_equidistant_frames
    )

    assert isinstance(fig, matplotlib.figure.Figure)
    if isinstance(superfeature_names, str):
        assert isinstance(axes, matplotlib.axes.Subplot)
    else:
        for ax in axes:
            assert isinstance(ax, matplotlib.axes.Subplot)


@pytest.mark.parametrize("superfeature_name", ["xxx", ["AR[4605,4607,4603,4606,4604]", "xxx"]])
def test_envpartners_occurrences_raises(dynophore, superfeature_name):

    with pytest.raises(KeyError):
        plot.static.envpartners_occurrences(dynophore, superfeature_name)


@pytest.mark.parametrize(
    "superfeature_names, kind",
    [
        ("AR[4605,4607,4603,4606,4604]", "line"),
        (["AR[4605,4607,4603,4606,4604]", "AR[4622,4615,4623,4613,4614,4621]"], "line"),
        ("AR[4605,4607,4603,4606,4604]", "hist"),
    ],
)
def test_envpartners_distances(dynophore, superfeature_names, kind):

    fig, axes = plot.static.envpartners_distances(dynophore, superfeature_names, kind)

    assert isinstance(fig, matplotlib.figure.Figure)
    if isinstance(superfeature_names, str):
        assert isinstance(axes, matplotlib.axes.Subplot)
    else:
        for ax in axes:
            assert isinstance(ax, matplotlib.axes.Subplot)


@pytest.mark.parametrize(
    "superfeature_name, kind", [("xxx", "line"), ("AR[4605,4607,4603,4606,4604]", "xxx")]
)
def test_envpartner_distances_raises(dynophore, superfeature_name, kind):

    with pytest.raises(KeyError):
        plot.static.envpartners_distances(dynophore, superfeature_name, kind)


@pytest.mark.parametrize(
    "superfeature_name, n_equidistant_frames",
    [
        ("AR[4605,4607,4603,4606,4604]", 1000),
    ],
)
def test_envpartners_all_in_one(dynophore, superfeature_name, n_equidistant_frames):

    fig, axes = plot.static.envpartners_all_in_one(
        dynophore, superfeature_name, n_equidistant_frames
    )
    print(axes.size)

    assert isinstance(fig, matplotlib.figure.Figure)
    assert axes.size == 4
    assert isinstance(axes[0][0], matplotlib.axes.Subplot)


@pytest.mark.parametrize("superfeature_name", ["xxx"])
def test_envpartners_all_in_one_raises(dynophore, superfeature_name):

    with pytest.raises(KeyError):
        plot.static.envpartners_all_in_one(dynophore, superfeature_name)