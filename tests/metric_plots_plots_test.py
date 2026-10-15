"""Direct tests for the low-level QC plotting functions."""

import matplotlib.pyplot as plt
import numpy as np
import polars as pl
import pytest
from plotnine import ggplot

from tsoppy.metric_plots.plots import (
    TABLEAU_20,
    Plot_bar_metric,
    Plot_contamination_scatter,
)


def _bar_data() -> pl.DataFrame:
    """Return representative bar-plot input."""
    return pl.DataFrame(
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


def _contamination_data() -> pl.DataFrame:
    """Return representative contamination-plot input."""
    return pl.DataFrame(
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
    data: pl.DataFrame | None = None,
    **kwargs,
):
    """Create a standard test bar plot."""
    if data is None:
        data = _bar_data()

    return Plot_bar_metric(
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
    data: pl.DataFrame | None = None,
    **kwargs,
):
    """Create a standard test contamination plot."""
    if data is None:
        data = _contamination_data()

    return Plot_contamination_scatter(
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
# Plot_bar_metric
# ---------------------------------------------------------------------------


def test_plot_bar_metric_draws():
    """A representative bar plot renders successfully."""
    plot = _basic_bar_plot()

    assert isinstance(plot, ggplot)

    figure = plot.draw()

    assert figure is not None
    assert len(figure.axes) >= 1

    plt.close(figure)


@pytest.mark.parametrize(
    "scale_name, column_name",
    [
        ("x", "PLOT_SAMPLE_ID"),
        ("fill", "PLOT_RUN"),
    ],
)
def test_plot_bar_metric_preserves_category_order(scale_name, column_name):
    """Sample and run legend categories retain their input/first-occurrence order."""
    data = _bar_data()

    plot = _basic_bar_plot(data)

    scale = plot.scales.get_scales(scale_name)

    assert list(scale.limits) == (
        data.get_column(column_name).unique(maintain_order=True).to_list()
    )


@pytest.mark.parametrize(
    ("guideline_kwargs", "want_layers"),
    [
        (
            # no guideline configured
            {},
            1,
        ),
        (
            # a configured guideline adds line and text annotation layers
            {
                "hline_y": 25.0,
                "hline_label": "USL_Guideline: 25",
            },
            3,
        ),
    ],
)
def test_plot_bar_metric_layer_count(guideline_kwargs, want_layers):
    """Plot layer count reflects whether a guideline is configured."""
    plot = _basic_bar_plot(**guideline_kwargs)

    assert len(plot.layers) == want_layers

    figure = plot.draw()
    assert figure is not None

    plt.close(figure)


@pytest.mark.parametrize(
    "guideline_kwargs",
    [
        (
            # explicit custom styling
            {
                "hline_y": 25.0,
                "hline_alpha": 0.25,
                "hline_color": "blue",
                "hline_size": 2.0,
                "hline_label": "Configured guideline",
                "ann_y_offset": 2.5,
            }
        ),
        (
            # None alpha/color values fall back to module guideline defaults
            {
                "hline_y": 25.0,
                "hline_alpha": None,
                "hline_color": None,
                "hline_label": "Guideline",
            }
        ),
    ],
)
def test_plot_bar_metric_draws_with_guideline_style(guideline_kwargs):
    """Guideline rendering accepts both explicit and default styling parameters."""
    plot = _basic_bar_plot(**guideline_kwargs)

    figure = plot.draw()

    assert figure is not None

    plt.close(figure)


@pytest.mark.parametrize(
    ("cart_ylim", "y_tick_step", "want_breaks"),
    [
        (
            # positive tick spacing creates explicit y-axis breaks
            None,
            10,
            [0.0, 10.0, 20.0, 30.0],
        ),
        (
            # configured upper plotting limit controls generated ticks
            (0, 50),
            10,
            [0.0, 10.0, 20.0, 30.0, 40.0, 50.0],
        ),
        (
            # missing tick spacing adds no y scale
            None,
            None,
            None,
        ),
        (
            # zero tick spacing adds no y scale
            None,
            0,
            None,
        ),
    ],
)
def test_plot_bar_metric_y_breaks(cart_ylim, y_tick_step, want_breaks):
    """Y-axis breaks follow the configured tick step and cartesian limit."""
    kwargs = {"y_tick_step": y_tick_step}

    if cart_ylim is not None:
        kwargs["cart_ylim"] = cart_ylim

    plot = _basic_bar_plot(**kwargs)

    y_scale = plot.scales.get_scales("y")

    if want_breaks is None:
        assert y_scale is None
    else:
        assert list(y_scale.breaks) == want_breaks


def test_plot_bar_metric_nan_max_does_not_fail_tick_generation():
    """All-missing numeric data does not create invalid y-axis breaks."""
    data = pl.DataFrame(
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


@pytest.mark.parametrize(
    "sample_count",
    [
        1,  # a single sample
        8,  # multiple run categories
    ],
)
def test_plot_bar_metric_draws_for_sample_count(sample_count):
    """Bar plotting renders for both a single sample and many run categories."""
    data = pl.DataFrame(
        {
            "PLOT_SAMPLE_ID": [
                f"{index:03d} | SAMPLE_{index}" for index in range(1, sample_count + 1)
            ],
            "VALUE": [float(index * 10) for index in range(1, sample_count + 1)],
            "PLOT_RUN": [
                f"{index:03d} | RUN_{index}" for index in range(1, sample_count + 1)
            ],
        }
    )

    plot = _basic_bar_plot(data)

    figure = plot.draw()

    assert figure is not None

    plt.close(figure)


# ---------------------------------------------------------------------------
# Plot_contamination_scatter
# ---------------------------------------------------------------------------


def test_plot_contamination_scatter_draws():
    """Representative contamination data renders successfully."""
    plot = _basic_contamination_plot()

    assert isinstance(plot, ggplot)

    figure = plot.draw()

    assert figure is not None
    assert len(figure.axes) >= 1

    plt.close(figure)


@pytest.mark.parametrize(
    ("label_var", "want_layers"),
    [
        (
            None,
            6,
        ),
        (
            # sample-label rendering adds one geom_text layer
            "contamination_label",
            7,
        ),
    ],
)
def test_plot_contamination_scatter_layer_count(label_var, want_layers):
    """Plot layer count reflects whether sample labels are configured."""
    plot = Plot_contamination_scatter(
        data=_contamination_data(),
        color_var="RUN_LABEL",
        label_var=label_var,
        guide_title="Run",
        title="Contamination test",
        max_contamination_score=5000,
        usl_contamination_score=1457,
        usl_contamination_pval=0.05,
    )

    assert len(plot.layers) == want_layers

    figure = plot.draw()
    assert figure is not None

    plt.close(figure)


def _custom_colors_contamination_plot():
    """Build a plot using explicit run-color mappings."""
    return Plot_contamination_scatter(
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


def _large_score_range_contamination_plot():
    """Build a plot with contamination scores above the default 5000 floor."""
    data = _contamination_data().with_columns(
        pl.Series(
            "DNA_CONTAMINATION_SCORE",
            [
                100.0,
                750.0,
                7200.0,
            ],
        )
    )

    return Plot_contamination_scatter(
        data=data,
        color_var="RUN_LABEL",
        label_var=None,
        guide_title="Run",
        title="Large contamination score",
        max_contamination_score=7200,
        usl_contamination_score=1457,
        usl_contamination_pval=0.05,
    )


def _single_sample_contamination_plot():
    """Build a plot from a single-sample dataset."""
    data = pl.DataFrame(
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

    return Plot_contamination_scatter(
        data=data,
        color_var="RUN_LABEL",
        label_var=None,
        guide_title="Run",
        title="Single sample",
        max_contamination_score=5000,
        usl_contamination_score=1457,
        usl_contamination_pval=0.05,
    )


@pytest.mark.parametrize(
    "build_plot",
    [
        _custom_colors_contamination_plot,
        _large_score_range_contamination_plot,
        _single_sample_contamination_plot,
    ],
)
def test_plot_contamination_scatter_draws_under_varied_conditions(build_plot):
    """Contamination scatter renders under varied data and styling conditions."""
    figure = build_plot().draw()

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

    assert list(color_scale.limits) == (
        data.get_column("RUN_LABEL").unique(maintain_order=True).to_list()
    )


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
    plot = Plot_contamination_scatter(
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
