"""Direct tests for the low-level QC plotting functions."""

import matplotlib.pyplot as plt
import numpy as np
import polars as pl
import pytest
from pytest import mark
from plotnine import ggplot

from tsoppy.metric_plots.plots import (
    Plot_bar_metric,
    Plot_contamination_scatter,
)

test_bar_plot_data = pl.DataFrame(
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


test_scatter_plot_data = pl.DataFrame(
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


@mark.parametrize(
    "inputs, wants",
    [
        (
            (test_bar_plot_data, None, None, None, None, None, None, None, None),
            (
                ggplot,
                1,
                "x",
                "PLOT_SAMPLE_ID",
                ["001 | SAMPLE_A", "001 | SAMPLE_B", "002 | SAMPLE_C"],
            ),
        ),
        (
            (
                pl.DataFrame(
                    {
                        "PLOT_SAMPLE_ID": ["001 | SAMPLE_A"],
                        "VALUE": [10.0],
                        "PLOT_RUN": ["001 | RUN_A"],
                    }
                ),
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
            ),
            (ggplot, 1, "x", "PLOT_SAMPLE_ID", ["001 | SAMPLE_A"]),
        ),
        (
            (
                test_bar_plot_data,
                0.25,
                "blue",
                "USL_Guideline: 25",
                2.0,
                25.0,
                2.5,
                None,
                None,
            ),
            (
                ggplot,
                3,
                "fill",
                "PLOT_RUN",
                True,
            ),
        ),
        (
            (
                test_bar_plot_data,
                None,
                None,
                None,
                None,
                None,
                None,
                10,
                None,
            ),
            (
                ggplot,
                1,
                "y",
                "PLOT_SAMPLE_ID",
                [0.0, 10.0, 20.0, 30.0],
            ),
        ),
        (
            (
                test_bar_plot_data,
                None,
                None,
                None,
                None,
                None,
                None,
                10,
                (0, 50),
            ),
            (
                ggplot,
                1,
                "y",
                "PLOT_SAMPLE_ID",
                [0.0, 10.0, 20.0, 30.0, 40.0, 50.0],
            ),
        ),
        (
            (
                test_bar_plot_data,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
            ),
            (
                ggplot,
                1,
                "y",
                "PLOT_SAMPLE_ID",
                None,
            ),
        ),
        (
            (
                test_bar_plot_data,
                None,
                None,
                None,
                None,
                None,
                None,
                0,
                None,
            ),
            (
                ggplot,
                1,
                "y",
                "PLOT_SAMPLE_ID",
                None,
            ),
        ),
        (
            (
                pl.DataFrame(
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
                ),
                None,
                None,
                None,
                None,
                None,
                None,
                10,
                None,
            ),
            (
                ggplot,
                1,
                "y",
                "PLOT_SAMPLE_ID",
                None,
            ),
        ),
    ],
)
def test_plot_bar_metric(inputs, wants):
    got = Plot_bar_metric(
        data=inputs[0],
        x_var="PLOT_SAMPLE_ID",
        y_var="VALUE",
        fill_var="PLOT_RUN",
        guide_title="Run",
        x_lab="Run index | Sample ID",
        y_lab="Value",
        title="Test metric",
        hline_alpha=inputs[1],
        hline_color=inputs[2],
        hline_label=inputs[3],
        hline_size=inputs[4],
        hline_y=inputs[5],
        ann_y_offset=inputs[6],
        y_tick_step=inputs[7],
        cart_ylim=inputs[8],
    )
    assert isinstance(got, wants[0])
    assert len(got.layers) == wants[1]
    scales = got.scales.get_scales(wants[2])
    if scales:
        if scales.limits:
            assert list(scales.limits) == (
                inputs[0].get_column(wants[3]).unique(maintain_order=True).to_list()
            )
        if scales.breaks:
            assert scales.breaks == wants[4]

    figure = got.draw()

    assert figure is not None
    assert len(figure.axes) >= 1

    plt.close(figure)


@pytest.mark.parametrize(
    "inputs, wants",
    [
        (
            (
                test_scatter_plot_data,
                "RUN_LABEL",
                None,
                "Run",
                "Contamination test",
                5000,
                1457,
                0.05,
                None,
            ),
            (
                ggplot,
                6,
                ["001 | RUN_A", "002 | RUN_B"],
            ),
        ),
        (
            (
                pl.DataFrame(
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
                ),
                "RUN_LABEL",
                None,
                "Run",
                "Single sample",
                5000,
                1000,
                0.1,
                None,
            ),
            (
                ggplot,
                6,
                ["001 | RUN_A"],
            ),
        ),
        (
            (
                test_scatter_plot_data,
                "RUN_LABEL",
                "contamination_label",
                "Run",
                "Contamination test",
                5000,
                5000,
                0.01,
                None,
            ),
            (
                ggplot,
                7,
                ["001 | RUN_A", "002 | RUN_B"],
            ),
        ),
        (
            (
                test_scatter_plot_data,
                "highlighted_run",
                None,
                "Latest run",
                "Custom colors",
                5000,
                1457,
                0.05,
                {
                    "False": "#888686",
                    "True": "#C02F2F",
                },
            ),
            (
                ggplot,
                6,
                ["True", "False"],
            ),
        ),
        (
            (
                test_scatter_plot_data.with_columns(
                    pl.Series(
                        "DNA_CONTAMINATION_SCORE",
                        [
                            100.0,
                            750.0,
                            7200.0,
                        ],
                    )
                ),
                "RUN_LABEL",
                None,
                "Run",
                "Large contamination score",
                7200,
                1457,
                0.05,
                None,
            ),
            (
                ggplot,
                6,
                ["001 | RUN_A", "002 | RUN_B"],
            ),
        ),
    ],
)
def test_plot_contamination_scatter(inputs, wants):
    got = Plot_contamination_scatter(
        data=inputs[0],
        color_var=inputs[1],
        label_var=inputs[2],
        guide_title=inputs[3],
        title=inputs[4],
        max_contamination_score=inputs[5],
        usl_contamination_score=inputs[6],
        usl_contamination_pval=inputs[7],
        color_values=inputs[8],
    )
    assert isinstance(got, wants[0])
    assert len(got.layers) == wants[1]
    scale_color = got.scales.get_scales("color")

    assert scale_color is not None
    assert list(scale_color.limits) == wants[2]

    figure = got.draw()

    assert figure is not None
    assert len(figure.axes) >= 1

    plt.close(figure)
