"""
Unit tests for dynophore.viz.plot.static.

Will only test if static plotting raises errors.
"""

import pytest
import matplotlib

from dynophores.viz import plot


@pytest.mark.parametrize(
    "superfeature_ids",
    [
        ("all"),  # Default
        (("all",)),
        ("AR[4605,4607,4603,4606,4604]"),
        (["AR[4605,4607,4603,4606,4604]", "AR[4622,4615,4623,4613,4614,4621]"]),
        (["all", "AR[4622,4615,4623,4613,4614,4621]"]),
    ],
)
def test_superfeatures_vs_envpartners(dynophore, superfeature_ids):

    fig, ax = plot.static.superfeatures_vs_envpartners(dynophore, superfeature_ids)

    assert isinstance(fig, matplotlib.figure.Figure)
    assert isinstance(ax, matplotlib.axes.Subplot)


@pytest.mark.parametrize("superfeature_ids", ["xxx"])
def test_superfeatures_vs_envpartners_raises(dynophore, superfeature_ids):

    with pytest.raises(KeyError):
        plot.static.superfeatures_vs_envpartners(dynophore, superfeature_ids)


@pytest.mark.parametrize(
    "superfeature_ids, color_by_feature_type, frames_range, frames_step_size",
    [
        ("all", True, [0, None], 1),  # Defaults
        ("all", False, [10, 100], 1),
        (("all",), False, [0, None], 100),
        ("AR[4605,4607,4603,4606,4604]", False, [0, None], 1),
        (
            ["AR[4605,4607,4603,4606,4604]", "AR[4622,4615,4623,4613,4614,4621]"],
            False,
            [0, None],
            1,
        ),
        (["all", "AR[4622,4615,4623,4613,4614,4621]"], False, [0, None], 1),
    ],
)
def test_superfeatures_occurrences(
    dynophore, superfeature_ids, color_by_feature_type, frames_range, frames_step_size
):

    fig, ax = plot.static.superfeatures_occurrences(
        dynophore, superfeature_ids, color_by_feature_type, frames_range, frames_step_size
    )

    assert isinstance(fig, matplotlib.figure.Figure)
    assert isinstance(ax, matplotlib.axes.Subplot)


@pytest.mark.parametrize("superfeature_ids", ["xxx"])
def test_superfeatures_occurrences_raises(dynophore, superfeature_ids):

    with pytest.raises(KeyError):
        plot.static.superfeatures_occurrences(dynophore, superfeature_ids)


@pytest.mark.parametrize(
    "superfeature_ids, frames_range, frames_step_size",
    [
        ("AR[4605,4607,4603,4606,4604]", [0, None], 1),
        (["AR[4605,4607,4603,4606,4604]", "AR[4622,4615,4623,4613,4614,4621]"], [0, None], 10),
        (["AR[4605,4607,4603,4606,4604]", "AR[4622,4615,4623,4613,4614,4621]"], [10, 90], 1),
        (["AR[4605,4607,4603,4606,4604]", "AR[4622,4615,4623,4613,4614,4621]"], [10, 90], 10),
    ],
)
def test_envpartners_occurrences(dynophore, superfeature_ids, frames_range, frames_step_size):

    fig, axes = plot.static.envpartners_occurrences(
        dynophore, superfeature_ids, frames_range, frames_step_size
    )

    assert isinstance(fig, matplotlib.figure.Figure)
    if isinstance(superfeature_ids, str):
        assert isinstance(axes, matplotlib.axes.Subplot)
    else:
        for ax in axes:
            assert isinstance(ax, matplotlib.axes.Subplot)


@pytest.mark.parametrize("superfeature_id", ["xxx", ["AR[4605,4607,4603,4606,4604]", "xxx"]])
def test_envpartners_occurrences_raises(dynophore, superfeature_id):

    with pytest.raises(KeyError):
        plot.static.envpartners_occurrences(dynophore, superfeature_id)


@pytest.mark.parametrize(
    "superfeature_ids, kind",
    [
        ("AR[4605,4607,4603,4606,4604]", "line"),
        (["AR[4605,4607,4603,4606,4604]", "AR[4622,4615,4623,4613,4614,4621]"], "line"),
        ("AR[4605,4607,4603,4606,4604]", "hist"),
    ],
)
def test_envpartners_distances(dynophore, superfeature_ids, kind):

    fig, axes = plot.static.envpartners_distances(dynophore, superfeature_ids, kind)

    assert isinstance(fig, matplotlib.figure.Figure)
    if isinstance(superfeature_ids, str):
        assert isinstance(axes, matplotlib.axes.Subplot)
    else:
        for ax in axes:
            assert isinstance(ax, matplotlib.axes.Subplot)


@pytest.mark.parametrize(
    "superfeature_id, kind", [("xxx", "line"), ("AR[4605,4607,4603,4606,4604]", "xxx")]
)
def test_envpartner_distances_raises(dynophore, superfeature_id, kind):

    with pytest.raises(KeyError):
        plot.static.envpartners_distances(dynophore, superfeature_id, kind)


@pytest.mark.parametrize(
    "superfeature_id, frames_range, frames_step_size",
    [
        ("AR[4605,4607,4603,4606,4604]", [0, None], 1),
        (("AR[4605,4607,4603,4606,4604]",), [0, None], 1),
    ],
)
def test_envpartners_all_in_one(dynophore, superfeature_id, frames_range, frames_step_size):

    fig, axes = plot.static.envpartners_all_in_one(
        dynophore, superfeature_id, frames_range, frames_step_size
    )
    print(axes.size)

    assert isinstance(fig, matplotlib.figure.Figure)
    assert axes.size == 4
    assert isinstance(axes[0][0], matplotlib.axes.Subplot)


@pytest.mark.parametrize("superfeature_id", ["xxx"])
def test_envpartners_all_in_one_raises(dynophore, superfeature_id):

    with pytest.raises(KeyError):
        plot.static.envpartners_all_in_one(dynophore, superfeature_id)
