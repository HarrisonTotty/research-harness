"""Tests for the blog-derived figure style."""

import matplotlib
import pytest
from matplotlib.colors import to_hex

from figures import style


def _relative_luminance(hex_color: str) -> float:
    def channel(value: int) -> float:
        c = value / 255
        return c / 12.92 if c <= 0.04045 else ((c + 0.055) / 1.055) ** 2.4

    r, g, b = (int(hex_color[i : i + 2], 16) for i in (1, 3, 5))
    return 0.2126 * channel(r) + 0.7152 * channel(g) + 0.0722 * channel(b)


def _contrast_on_paper(hex_color: str) -> float:
    return 1.05 / (_relative_luminance(hex_color) + 0.05)


def test_categorical_slots_are_six_unique_hexes():
    assert len(style.CATEGORICAL) == 6
    assert len(set(style.CATEGORICAL)) == 6
    for color in style.CATEGORICAL:
        assert color.startswith("#")
        assert len(color) == 7


@pytest.mark.parametrize("color", style.CATEGORICAL)
def test_categorical_slot_clears_3_to_1_contrast_on_paper(color):
    assert _contrast_on_paper(color) >= 3.0


def test_sequential_anchors_darken_monotonically():
    luminances = [_relative_luminance(c) for c in style.SEQUENTIAL_ANCHORS]
    assert luminances == sorted(luminances, reverse=True)


def test_ordinal_steps_darken_monotonically_from_a_readable_light_end():
    luminances = [_relative_luminance(c) for c in style.ORDINAL_STEPS]
    assert luminances == sorted(luminances, reverse=True)
    assert _contrast_on_paper(style.ORDINAL_STEPS[0]) >= 2.0


def test_diverging_anchors_peak_in_lightness_at_the_midpoint():
    luminances = [_relative_luminance(c) for c in style.DIVERGING_ANCHORS]
    mid = len(luminances) // 2
    assert luminances[mid] == max(luminances)
    assert luminances[:mid] == sorted(luminances[:mid])
    assert luminances[mid:] == sorted(luminances[mid:], reverse=True)


def test_rc_params_only_names_known_rcparams():
    unknown = set(style.rc_params()) - set(matplotlib.rcParams)
    assert unknown == set()


def test_context_applies_and_restores_rcparams():
    before = matplotlib.rcParams["grid.color"]

    with style.context():
        assert matplotlib.rcParams["grid.color"] == style.PARCHMENT
        assert matplotlib.rcParams["axes.edgecolor"] == style.MIST

    assert matplotlib.rcParams["grid.color"] == before


def test_context_installs_the_categorical_property_cycle():
    with style.context():
        cycle = matplotlib.rcParams["axes.prop_cycle"].by_key()["color"]

    assert tuple(cycle) == style.CATEGORICAL


def test_register_colormaps_is_idempotent():
    style.register_colormaps()
    style.register_colormaps()

    assert "blog-sequential" in matplotlib.colormaps
    assert "blog-diverging" in matplotlib.colormaps
    assert "blog-categorical" in matplotlib.colormaps


def test_categorical_cmap_lists_the_slots_in_order():
    sampled = [
        to_hex(style.CATEGORICAL_CMAP(i)) for i in range(style.CATEGORICAL_CMAP.N)
    ]
    assert sampled == list(style.CATEGORICAL)
