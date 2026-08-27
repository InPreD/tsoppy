"""Direct tests for the low-level QC plotting functions."""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from plotnine import ggplot

from tsoppy.metric_plots.plots import (
    HV_LINE_ALPHA,
    HV_LINE_COLOR,
    TABLEAU_20,
    plot_contamination_scatter,
    plot_tsoppy_barplot,
)


def _bar_data() -> pd.DataFrame:
    """Return representative bar-plot input."""
    return pd.DataFrame(
        {
            "PLOT_SAMPLE_ID": [
                "001 | SAMPLE_A",
                "001 | SAMPLE_B",
                "002 | SAMPLE_C",
            ],
            "VALUE": [
                10.0,
                20.0,
                30.0,
            ],
            "PLOT_RUN": [
                "001 | RUN_A",
                "001 | RUN_A",
                "002 | RUN_B",
            ],
        }
    )


def _contamination_data() -> pd.DataFrame:
    """Return representative contamination-plot input."""
    return pd.DataFrame(
        {
            "DNA_CONTAMINATION_SCORE": [
                100.0,
                750.0,
                2000.0,
            ],
            "DNA_CONTAMINATION_P_VALUE": [
                0.01,
                0.08,
                0.20,
            ],
            "RUN_LABEL": [
                "001 | RUN_A",
                "001 | RUN_A",
                "002 | RUN_B",
            ],
            "highlighted_run": [
                "True",
                "True",
                "False",
            ],
            "contamination_label": [
                "SAMPLE_A",
                "SAMPLE_B",
                "",
            ],
        }
    )


def _basic_bar_plot(
    data: pd.DataFrame | None = None,
    **kwargs,
):
    """Create a standard test bar plot."""
    if data is None:
        data = _bar_data()

    return plot_tsoppy_barplot(
        data=data,
        x_var="PLOT_SAMPLE_ID",
        y_var="VALUE",
        fill_var="PLOT_RUN",
        guide_title="Run",
        x_lab="Run index | Sample ID",
        y_lab="Value",
        title="Test metric",
        **kwargs,
    )


def _basic_contamination_plot(
    data: pd.DataFrame | None = None,
    **kwargs,
):
    """Create a standard test contamination plot."""
    if data is None:
        data = _contamination_data()

    return plot_contamination_scatter(
        data=data,
        color_var="RUN_LABEL",
        label_var=None,
        guide_title="Run",
        title="Contamination test",
        max_contamination_score=5000,
        usl_contamination_score=1457,
        usl_contamination_pval=0.05,
        **kwargs,
    )


# ---------------------------------------------------------------------------
# plot_tsoppy_barplot
# ---------------------------------------------------------------------------


def test_plot_tsoppy_barplot_returns_ggplot():
    """Bar plotting helper returns a plotnine ggplot object."""
    plot = _basic_bar_plot()

    assert isinstance(plot, ggplot)


def test_plot_tsoppy_barplot_draws():
    """A representative bar plot renders successfully."""
    plot = _basic_bar_plot()

    figure = plot.draw()

    assert figure is not None
    assert len(figure.axes) >= 1

    plt.close(figure)


def test_plot_tsoppy_barplot_preserves_x_category_order():
    """Sample categories retain their input ordering."""
    data = _bar_data()

    plot = _basic_bar_plot(data)

    x_scale = plot.scales.get_scales("x")

    assert list(x_scale.limits) == data["PLOT_SAMPLE_ID"].drop_duplicates().tolist()


def test_plot_tsoppy_barplot_preserves_fill_category_order():
    """Run legend categories retain their first-occurrence order."""
    data = _bar_data()

    plot = _basic_bar_plot(data)

    fill_scale = plot.scales.get_scales("fill")

    assert list(fill_scale.limits) == data["PLOT_RUN"].drop_duplicates().tolist()


def test_plot_tsoppy_barplot_without_guideline_has_one_layer():
    """Basic bar plot contains only the bar layer."""
    plot = _basic_bar_plot()

    assert len(plot.layers) == 1


def test_plot_tsoppy_barplot_adds_guideline_layers():
    """Horizontal guideline adds line and text annotation layers."""
    plot = _basic_bar_plot(
        hline_y=25.0,
        hline_label="USL_Guideline: 25",
    )

    assert len(plot.layers) == 3

    figure = plot.draw()
    assert figure is not None

    plt.close(figure)


def test_plot_tsoppy_barplot_guideline_accepts_custom_style():
    """Guideline rendering accepts configured styling parameters."""
    plot = _basic_bar_plot(
        hline_y=25.0,
        hline_alpha=0.25,
        hline_color="blue",
        hline_size=2.0,
        hline_label="Configured guideline",
        ann_y_offset=2.5,
    )

    figure = plot.draw()

    assert figure is not None

    plt.close(figure)


def test_plot_tsoppy_barplot_guideline_none_style_uses_defaults():
    """None alpha/color values fall back to module guideline defaults."""
    plot = _basic_bar_plot(
        hline_y=25.0,
        hline_alpha=None,
        hline_color=None,
        hline_label="Guideline",
    )

    figure = plot.draw()

    assert figure is not None
    assert HV_LINE_ALPHA is not None
    assert HV_LINE_COLOR is not None

    plt.close(figure)


def test_plot_tsoppy_barplot_adds_requested_y_breaks():
    """Positive y_tick_step creates explicit y-axis breaks."""
    plot = _basic_bar_plot(
        y_tick_step=10,
    )

    y_scale = plot.scales.get_scales("y")

    assert list(y_scale.breaks) == [
        0.0,
        10.0,
        20.0,
        30.0,
    ]


def test_plot_tsoppy_barplot_y_breaks_use_cartesian_upper_limit():
    """Configured upper plotting limit controls generated ticks."""
    plot = _basic_bar_plot(
        cart_ylim=(0, 50),
        y_tick_step=10,
    )

    y_scale = plot.scales.get_scales("y")

    assert list(y_scale.breaks) == [
        0.0,
        10.0,
        20.0,
        30.0,
        40.0,
        50.0,
    ]


@pytest.mark.parametrize(
    "tick_step",
    [
        None,
        0,
        -1,
    ],
)
def test_plot_tsoppy_barplot_non_positive_tick_step_adds_no_y_scale(
    tick_step,
):
    """Missing or non-positive tick spacing does not add a y scale."""
    plot = _basic_bar_plot(
        y_tick_step=tick_step,
    )

    assert plot.scales.get_scales("y") is None


def test_plot_tsoppy_barplot_nan_max_does_not_fail_tick_generation():
    """All-missing numeric data does not create invalid y-axis breaks."""
    data = pd.DataFrame(
        {
            "PLOT_SAMPLE_ID": [
                "001 | SAMPLE_A",
            ],
            "VALUE": [
                np.nan,
            ],
            "PLOT_RUN": [
                "001 | RUN_A",
            ],
        }
    )

    plot = _basic_bar_plot(
        data,
        y_tick_step=10,
    )

    assert isinstance(plot, ggplot)


def test_plot_tsoppy_barplot_single_sample_draws():
    """Bar plotting also works with only one sample."""
    data = pd.DataFrame(
        {
            "PLOT_SAMPLE_ID": [
                "001 | SAMPLE_A",
            ],
            "VALUE": [
                10.0,
            ],
            "PLOT_RUN": [
                "001 | RUN_A",
            ],
        }
    )

    plot = _basic_bar_plot(data)

    figure = plot.draw()

    assert figure is not None

    plt.close(figure)


def test_plot_tsoppy_barplot_many_runs_draws():
    """Multiple run categories can be rendered."""
    data = pd.DataFrame(
        {
            "PLOT_SAMPLE_ID": [
                f"{index:03d} | SAMPLE_{index}" for index in range(1, 9)
            ],
            "VALUE": [float(index * 10) for index in range(1, 9)],
            "PLOT_RUN": [f"{index:03d} | RUN_{index}" for index in range(1, 9)],
        }
    )

    plot = _basic_bar_plot(data)

    figure = plot.draw()

    assert figure is not None

    plt.close(figure)


# ---------------------------------------------------------------------------
# plot_contamination_scatter
# ---------------------------------------------------------------------------


def test_plot_contamination_scatter_returns_ggplot():
    """Contamination helper returns a plotnine ggplot object."""
    plot = _basic_contamination_plot()

    assert isinstance(plot, ggplot)


def test_plot_contamination_scatter_draws():
    """Representative contamination data renders successfully."""
    plot = _basic_contamination_plot()

    figure = plot.draw()

    assert figure is not None
    assert len(figure.axes) >= 1

    plt.close(figure)


def test_plot_contamination_scatter_without_labels_has_six_layers():
    """Base contamination plot contains expected geometry layers."""
    plot = _basic_contamination_plot()

    assert len(plot.layers) == 6


def test_plot_contamination_scatter_with_labels_adds_layer():
    """Sample-label rendering adds one geom_text layer."""
    plot = plot_contamination_scatter(
        data=_contamination_data(),
        color_var="highlighted_run",
        label_var="contamination_label",
        guide_title="Latest run",
        title="Highlighted run",
        max_contamination_score=5000,
        usl_contamination_score=1457,
        usl_contamination_pval=0.05,
        color_values={
            "False": "#888686",
            "True": "#C02F2F",
        },
    )

    assert len(plot.layers) == 7

    figure = plot.draw()

    assert figure is not None

    plt.close(figure)


def test_plot_contamination_scatter_custom_colors_draw():
    """Explicit run-color mappings are supported."""
    plot = plot_contamination_scatter(
        data=_contamination_data(),
        color_var="highlighted_run",
        label_var=None,
        guide_title="Latest run",
        title="Custom colors",
        max_contamination_score=5000,
        usl_contamination_score=1457,
        usl_contamination_pval=0.05,
        color_values={
            "False": "#888686",
            "True": "#C02F2F",
        },
    )

    figure = plot.draw()

    assert figure is not None

    plt.close(figure)


def test_plot_contamination_scatter_uses_default_palette():
    """Default contamination plot uses the shared Tableau palette."""
    plot = _basic_contamination_plot()

    color_scale = plot.scales.get_scales("color")

    assert color_scale is not None
    assert len(TABLEAU_20) > 0


def test_plot_contamination_scatter_preserves_run_order():
    """Run legend categories retain first-occurrence ordering."""
    data = _contamination_data()

    plot = _basic_contamination_plot(data)

    color_scale = plot.scales.get_scales("color")

    assert list(color_scale.limits) == data["RUN_LABEL"].drop_duplicates().tolist()


def test_plot_contamination_scatter_large_score_range_draws():
    """Contamination scores above 5000 can be plotted."""
    data = _contamination_data().copy()

    data.loc[
        2,
        "DNA_CONTAMINATION_SCORE",
    ] = 7200.0

    plot = plot_contamination_scatter(
        data=data,
        color_var="RUN_LABEL",
        label_var=None,
        guide_title="Run",
        title="Large contamination score",
        max_contamination_score=7200,
        usl_contamination_score=1457,
        usl_contamination_pval=0.05,
    )

    figure = plot.draw()

    assert figure is not None

    plt.close(figure)


def test_plot_contamination_scatter_single_sample_draws():
    """Contamination scatter supports a single sample."""
    data = pd.DataFrame(
        {
            "DNA_CONTAMINATION_SCORE": [
                100.0,
            ],
            "DNA_CONTAMINATION_P_VALUE": [
                0.01,
            ],
            "RUN_LABEL": [
                "001 | RUN_A",
            ],
        }
    )

    plot = plot_contamination_scatter(
        data=data,
        color_var="RUN_LABEL",
        label_var=None,
        guide_title="Run",
        title="Single sample",
        max_contamination_score=5000,
        usl_contamination_score=1457,
        usl_contamination_pval=0.05,
    )

    figure = plot.draw()

    assert figure is not None

    plt.close(figure)


@pytest.mark.parametrize(
    "score_limit,pvalue_limit",
    [
        (1457.0, 0.05),
        (1000.0, 0.10),
        (5000.0, 0.01),
    ],
)
def test_plot_contamination_scatter_guideline_values_draw(
    score_limit,
    pvalue_limit,
):
    """Different contamination guideline combinations render."""
    plot = plot_contamination_scatter(
        data=_contamination_data(),
        color_var="RUN_LABEL",
        label_var=None,
        guide_title="Run",
        title="Guideline test",
        max_contamination_score=5000,
        usl_contamination_score=score_limit,
        usl_contamination_pval=pvalue_limit,
    )

    figure = plot.draw()

    assert figure is not None

    plt.close(figure)
