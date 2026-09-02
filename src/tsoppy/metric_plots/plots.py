"""Low-level Plotnine helpers for metric QC plots."""

import numpy as np
import polars as pl
from plotnine import (
    aes,
    annotate,
    coord_cartesian,
    element_line,
    element_text,
    geom_col,
    geom_hline,
    geom_point,
    geom_text,
    geom_vline,
    ggplot,
    ggtitle,
    guide_legend,
    guides,
    scale_color_manual,
    scale_fill_manual,
    scale_x_continuous,
    scale_x_discrete,
    scale_y_continuous,
    theme,
    theme_minimal,
    xlab,
    ylab,
)

# Shared styling for plot axes and lines.
AXIS_LINE_COLOR = "#888686"

# Styling for the horizontal/vertical guideline markers (LSL/USL) drawn on plots.
HV_LINE_COLOR = "#e00202"
HV_LINE_SIZE = 1
HV_LINE_ALPHA = 0.75

# Categorical palette (Tableau 20, first 10 colors) used to fill/color
# series such as runs or samples across all QC plots.
TABLEAU_20 = [
    "#4E79A7",
    "#F28E2B",
    "#E15759",
    "#76B7B2",
    "#59A14F",
    "#EDC948",
    "#B07AA1",
    "#FF9DA7",
    "#9C755F",
    "#BAB0AC",
]


def Plot_bar_metric(
    data: pl.DataFrame | None = None,
    x_var: str | None = None,
    y_var: str | None = None,
    fill_var: str | None = None,
    guide_title: str | None = None,
    x_lab: str | None = None,
    y_lab: str | None = None,
    cart_ylim: tuple | None = None,
    title: str | None = None,
    x_lab_angle: int = 45,
    fig_size: tuple = (15, 6),
    alpha_value: float = 0.8,
    hline_y: float | None = None,
    hline_alpha: float = HV_LINE_ALPHA,
    hline_color: str = HV_LINE_COLOR,
    hline_label: str = "",
    hline_size: float = HV_LINE_SIZE,
    ann_y_offset: float = 0.0,
    y_tick_step: int | float | None = None,
) -> ggplot:
    """Create a reusable QC metric bar plot."""

    plot = (
        ggplot(data, aes(x=x_var, y=y_var))
        + geom_col(aes(fill=fill_var), alpha=alpha_value)
        + scale_fill_manual(
            values=TABLEAU_20,
            limits=(data.get_column(fill_var).unique(maintain_order=True).to_list()),
        )
        + guides(fill=guide_legend(guide_title))
        + xlab(x_lab)
        + ylab(y_lab)
        + coord_cartesian(ylim=cart_ylim)
        + ggtitle(title)
        + theme_minimal()
        + theme(
            axis_text_x=element_text(
                angle=x_lab_angle,
                size=7,
                ha="right",
                hjust=1.5,
                family="monospace",
            ),
            axis_ticks_major_x=element_line(
                color=AXIS_LINE_COLOR,
            ),
            figure_size=fig_size,
            axis_line=element_line(
                color=AXIS_LINE_COLOR,
            ),
            plot_title=element_text(
                weight="bold",
                ha="left",
            ),
        )
        + scale_x_discrete(
            limits=(data.get_column(x_var).unique(maintain_order=True).to_list()),
            breaks=(data.get_column(x_var).unique(maintain_order=True).to_list()),
            expand=(0.05, 0, 0.1, 0),
        )
    )

    if y_tick_step is not None and y_tick_step > 0:
        data_max = data.select(
            pl.col(y_var).cast(pl.Float64, strict=False).max()
        ).item()

        upper_lim = (
            cart_ylim[1]
            if (cart_ylim is not None and cart_ylim[1] is not None)
            else data_max
        )

        if upper_lim is not None and np.isfinite(upper_lim):
            upper_tick = upper_lim + y_tick_step

            y_breaks = list(
                np.arange(
                    0,
                    upper_tick,
                    y_tick_step,
                ).round(decimals=10)
            )

            plot = plot + scale_y_continuous(
                breaks=y_breaks,
            )

        if upper_lim is not None and np.isfinite(upper_lim):
            upper_tick = upper_lim + y_tick_step

            y_breaks = list(
                np.arange(
                    0,
                    upper_tick,
                    y_tick_step,
                ).round(decimals=10)
            )

            plot = plot + scale_y_continuous(
                breaks=y_breaks,
            )

    if hline_y is not None:
        plot = (
            plot
            + geom_hline(
                yintercept=hline_y,
                alpha=(HV_LINE_ALPHA if hline_alpha is None else hline_alpha),
                color=(HV_LINE_COLOR if hline_color is None else hline_color),
                size=hline_size,
            )
            + annotate(
                "text",
                len(data) + 1.25,
                hline_y + ann_y_offset,
                label=hline_label,
                angle=90,
                size=12,
            )
        )

    return plot


def Plot_contamination_scatter(
    data: pl.DataFrame | None = None,
    color_var: str | None = None,
    label_var: str | None = None,
    guide_title: str = "Run",
    title: str | None = None,
    x_lab_angle: int = 45,
    fig_size: tuple = (15, 6),
    max_contamination_score: float = 5000,
    usl_contamination_score: float | None = None,
    usl_contamination_pval: float | None = None,
    color_values: list | dict | None = None,
) -> ggplot:
    """Create a contamination assessment scatter plot."""

    x_breaks = [
        index * 1000 for index in range(int(max_contamination_score / 1000) + 2)
    ]

    y_breaks = [index * 0.1 for index in range(11)]

    selected_colors = color_values if color_values is not None else TABLEAU_20

    plot = (
        ggplot(
            data,
            aes(
                x="DNA_CONTAMINATION_SCORE",
                y="DNA_CONTAMINATION_P_VALUE",
            ),
        )
        + geom_point(
            aes(color=color_var),
            size=3,
        )
        + scale_color_manual(
            values=selected_colors,
            limits=(data.get_column(color_var).unique(maintain_order=True).to_list()),
        )
        + guides(color=guide_legend(guide_title))
        + xlab("Contamination score")
        + ylab("Contamination P-value")
        + coord_cartesian(
            xlim=(
                0,
                max_contamination_score + 500,
            ),
            ylim=(0.0, 1.05),
        )
        + scale_x_continuous(
            breaks=x_breaks,
            minor_breaks=[
                index * 500 for index in range(int(max_contamination_score / 500) + 2)
            ],
        )
        + scale_y_continuous(
            breaks=y_breaks,
            minor_breaks=[index * 0.05 for index in range(21)],
        )
        + ggtitle(title)
        + annotate(
            "rect",
            xmin=usl_contamination_score,
            xmax=max_contamination_score + 250,
            ymin=usl_contamination_pval,
            ymax=1,
            fill=HV_LINE_COLOR,
            alpha=0.3,
        )
        + geom_hline(
            yintercept=usl_contamination_pval,
            alpha=HV_LINE_ALPHA,
            color=HV_LINE_COLOR,
            size=HV_LINE_SIZE,
        )
        + annotate(
            "text",
            x=max_contamination_score + 500,
            y=usl_contamination_pval + 0.1,
            label=f"USL_Guideline: {usl_contamination_pval}",
            angle=90,
            size=10,
        )
        + geom_vline(
            xintercept=usl_contamination_score,
            alpha=HV_LINE_ALPHA,
            color=HV_LINE_COLOR,
            size=HV_LINE_SIZE,
        )
        + annotate(
            "text",
            x=usl_contamination_score + 200,
            y=1.05,
            label=f"USL_Guideline: {usl_contamination_score}",
            size=10,
        )
        + theme_minimal()
        + theme(
            axis_text_x=element_text(
                angle=x_lab_angle,
                vjust=0.5,
                hjust=1,
                size=7,
                family="monospace",
            ),
            figure_size=fig_size,
            axis_line=element_line(
                color=AXIS_LINE_COLOR,
            ),
            plot_title=element_text(
                weight="bold",
                ha="left",
            ),
            panel_grid_major=element_line(
                color=AXIS_LINE_COLOR,
                size=0.75,
                alpha=0.5,
            ),
            panel_grid_minor=element_line(
                color=AXIS_LINE_COLOR,
                size=0.25,
                alpha=0.5,
            ),
        )
    )

    if label_var is not None:
        plot = plot + geom_text(
            aes(label=label_var),
            angle=-25,
            nudge_y=0.02,
            size=7,
        )

    return plot
