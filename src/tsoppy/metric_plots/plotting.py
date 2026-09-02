"""Prepare and render workflow-specific metric QC plots."""

import logging
from pathlib import Path

import matplotlib.pyplot as plt
import polars as pl
from matplotlib.backends.backend_pdf import PdfPages
from plotnine import (
    aes,
    annotate,
    coord_cartesian,
    element_line,
    element_text,
    geom_hline,
    geom_point,
    geom_text,
    ggplot,
    ggtitle,
    guide_legend,
    guides,
    scale_color_manual,
    scale_x_continuous,
    scale_y_continuous,
    theme,
    theme_minimal,
    xlab,
    ylab,
)

from tsoppy.metric_plots.plots import (
    AXIS_LINE_COLOR,
    HV_LINE_ALPHA,
    HV_LINE_COLOR,
    HV_LINE_SIZE,
    TABLEAU_20,
    plot_bar_metric,
    plot_contamination_scatter,
)
from tsoppy.metric_plots.specs.plot_specs_workflows import PLOT_SPECS

logger = logging.getLogger(__name__)

ANGLE_X_NAMES = 90

_REQUIRED_COMMON_FIELDS = {
    "plot_kind",
    "source",
    "title",
}

_REQUIRED_FIELDS_BY_KIND = {
    "bar": {
        "x_var",
        "y_var",
        "fill_var",
        "x_lab",
        "y_lab",
        "value_spec",
    },
    "cluster_density_scatter": set(),
    "contamination_scatter": {
        "color_var",
        "label_var",
    },
}

SUPPORTED_WORKFLOWS = {
    "dragen",
    "localapp",
}

_REQUIRED_WORKFLOW_KEYS = {
    "plot",
    "index",
}


def resolve_plot_title(spec: dict, workflow: str) -> str:
    """Return the workflow-specific or shared plot title."""
    title = spec["title"]

    if isinstance(title, dict):
        return title[workflow]

    return title


def validate_plot_specs(plot_specs: dict) -> None:
    """Validate the structure of all plot specifications."""

    valid_kinds = set(_REQUIRED_FIELDS_BY_KIND)
    structural_errors = []

    for workflow in SUPPORTED_WORKFLOWS:
        seen_indexes = {}

        for spec_name, spec in plot_specs.items():
            workflow_entry = spec.get(workflow)

            if not isinstance(workflow_entry, dict):
                continue

            index = workflow_entry.get("index")

            if isinstance(index, int) and index >= 1:
                if index in seen_indexes:
                    message = (
                        f"Workflow '{workflow}' has duplicate plot index {index} "
                        f"in '{seen_indexes[index]}' and '{spec_name}'."
                    )
                    logger.error(message)
                    structural_errors.append(message)
                else:
                    seen_indexes[index] = spec_name

    if structural_errors:
        raise KeyError(
            f"Found {len(structural_errors)} duplicate plot index error(s):\n"
            + "\n".join(structural_errors)
        )

    for spec_name, spec in plot_specs.items():
        spec_errors = []

        for workflow in SUPPORTED_WORKFLOWS:
            if workflow not in spec:
                spec_errors.append(f"Missing workflow routing key: '{workflow}'.")
                continue

            workflow_entry = spec[workflow]

            missing_workflow_fields = _REQUIRED_WORKFLOW_KEYS - set(workflow_entry)

            if missing_workflow_fields:
                spec_errors.append(
                    f"Workflow '{workflow}' is missing fields: "
                    f"{', '.join(sorted(missing_workflow_fields))}."
                )
                continue

            if not isinstance(workflow_entry["plot"], bool):
                spec_errors.append(
                    f"Workflow '{workflow}': 'plot' must be bool, "
                    f"got {type(workflow_entry['plot']).__name__}."
                )

            if (
                not isinstance(workflow_entry["index"], int)
                or workflow_entry["index"] < 0
            ):
                spec_errors.append(
                    f"Workflow '{workflow}': 'index' must be a "
                    f"non-negative integer, got {workflow_entry['index']!r}."
                )

        required_fields = _REQUIRED_COMMON_FIELDS.copy()

        plot_kind = spec.get("plot_kind")

        if plot_kind in valid_kinds:
            required_fields.update(_REQUIRED_FIELDS_BY_KIND[plot_kind])
        elif plot_kind is not None:
            spec_errors.append(
                f"plot_kind '{plot_kind}' is not recognized. "
                f"Valid plot kinds: {', '.join(sorted(valid_kinds))}."
            )

        missing_fields = required_fields - set(spec)

        if missing_fields:
            spec_errors.append(f"Missing fields: {', '.join(sorted(missing_fields))}.")

        if spec_errors:
            message = f"Invalid plot specification '{spec_name}': " + " ".join(
                spec_errors
            )
            logger.error(message)
            structural_errors.append(message)

    if structural_errors:
        raise KeyError(
            f"Found {len(structural_errors)} plot specification error(s):\n"
            + "\n".join(structural_errors)
        )


def valid_metric_expr(column_name: str) -> pl.Expr:
    """Treat null values and literal NA strings as missing."""
    return pl.col(column_name).is_not_null() & (
        pl.col(column_name).cast(pl.Utf8) != "NA"
    )


def build_filter_expression(filter_spec: dict) -> pl.Expr:
    """Translate a plot filter specification into a Polars expression."""

    column_expr = pl.col(filter_spec["column"])

    if "contains" in filter_spec:
        return column_expr.cast(pl.Utf8).str.contains(filter_spec["contains"])

    if "equals" in filter_spec:
        return column_expr == filter_spec["equals"]

    if "not_equals" in filter_spec:
        return column_expr != filter_spec["not_equals"]

    raise ValueError(f"Unsupported filter specification: {filter_spec}")


def build_value_expression(
    value_spec: dict,
    alias_name: str | None = None,
) -> pl.Expr:
    """Create a raw, scaled, or derived plotting expression."""

    operation = value_spec.get("operation", "cast")

    if operation == "cast":
        expression = pl.col(value_spec["column"])

        if value_spec.get("dtype") is not None:
            expression = expression.cast(value_spec["dtype"])

    elif operation == "divide":
        expression = (
            pl.col(value_spec["column"]).cast(value_spec.get("dtype", pl.Float64))
            / value_spec["divisor"]
        )

    elif operation == "ratio":
        expression = pl.col(value_spec["numerator"]).cast(
            value_spec.get("dtype", pl.Float64)
        )

        numerator_divisor = value_spec.get("numerator_divisor")

        if numerator_divisor is not None:
            expression = expression / numerator_divisor

        expression = expression / pl.col(value_spec["denominator"]).cast(
            value_spec.get("dtype", pl.Float64)
        )

    else:
        raise ValueError(f"Unsupported value operation: {operation}")

    if alias_name is not None:
        expression = expression.alias(alias_name)

    return expression


def get_guideline_value(
    tables: dict,
    guideline_spec: dict,
) -> dict | None:
    """Extract one guideline value and its plotting parameters."""

    table = tables.get(guideline_spec["table"])

    if table is None or table.is_empty():
        return None

    id_column = guideline_spec.get(
        "id_column",
        "SAMPLE_ID",
    )

    if id_column not in table.columns:
        return None

    filtered_table = table.filter(pl.col(id_column) == guideline_spec["sample_id"])

    if filtered_table.is_empty():
        return None

    expression = build_value_expression(guideline_spec["value_spec"])

    try:
        value = filtered_table.select(expression).head(1).item()

    except (
        pl.exceptions.ColumnNotFoundError,
        pl.exceptions.InvalidOperationError,
        pl.exceptions.ComputeError,
    ):
        return None

    if value is None:
        return None

    if isinstance(value, str):
        normalized_value = value.strip().upper()

        if normalized_value in {
            "",
            "NA",
            "N/A",
            "NULL",
            "NONE",
        }:
            return None

    python_cast = guideline_spec.get("python_cast")

    if python_cast is not None:
        try:
            value = python_cast(value)
        except TypeError, ValueError:
            return None
    if guideline_spec.get("sample_id") == "LSL_Guideline" and float(value) == 0.0:
        return None

    return {
        "value": value,
        "label": guideline_spec.get(
            "label",
            f"{guideline_spec['label_prefix']}: {value}",
        ),
        "alpha": guideline_spec.get("alpha"),
        "color": guideline_spec.get("color"),
        "ann_y_offset": guideline_spec.get(
            "ann_y_offset",
            0,
        ),
    }


def get_available_guidelines(
    tables: dict,
    spec: dict,
) -> list[dict]:
    """Return each available LSL and USL guideline exactly once."""

    source_config = {
        "joint_qc_table": (
            "joint_qc_guideline_table",
            "RUN_ID",
        ),
        "dna_data_table": (
            "dna_guideline_table",
            "SAMPLE_ID",
        ),
        "rna_data_table": (
            "rna_guideline_table",
            "SAMPLE_ID",
        ),
        "data_table": (
            "guideline_table",
            "SAMPLE_ID",
        ),
        "merged_tables": (
            "guideline_table",
            "SAMPLE_ID",
        ),
    }

    explicit_specs = []

    explicit_guideline = spec.get("guideline")

    if isinstance(
        explicit_guideline,
        dict,
    ):
        explicit_specs.append(explicit_guideline)

    extra_guidelines = spec.get(
        "guidelines",
        [],
    )

    if extra_guidelines:
        explicit_specs.extend(extra_guidelines)

    explicit_by_sample_id = {
        guideline_spec["sample_id"]: guideline_spec
        for guideline_spec in explicit_specs
        if guideline_spec.get("sample_id") is not None
    }

    guidelines = []
    handled_sample_ids = set()

    auto_source = source_config.get(spec["source"])

    if auto_source is not None:
        table_name, id_column = auto_source

        if table_name in tables:
            for sample_id in (
                "LSL_Guideline",
                "USL_Guideline",
            ):
                guideline_spec = {
                    "table": table_name,
                    "id_column": id_column,
                    "sample_id": sample_id,
                    "value_spec": spec["value_spec"],
                    "python_cast": float,
                    "label_prefix": (sample_id),
                }

                explicit_override = explicit_by_sample_id.get(sample_id)

                if explicit_override is not None:
                    guideline_spec.update(explicit_override)

                guideline = get_guideline_value(
                    tables,
                    guideline_spec,
                )

                if guideline is None:
                    continue

                guideline["sample_id"] = sample_id

                guidelines.append(guideline)

                handled_sample_ids.add(sample_id)

    # Preserve explicitly configured guidelines
    # that were not already handled above.
    # This includes custom values such as
    # "Internal Guideline".
    for guideline_spec in explicit_specs:
        sample_id = guideline_spec.get("sample_id")

        if sample_id is None:
            continue

        if sample_id in handled_sample_ids:
            continue

        guideline = get_guideline_value(
            tables,
            guideline_spec,
        )

        if guideline is None:
            continue

        guideline["sample_id"] = sample_id

        guidelines.append(guideline)

        handled_sample_ids.add(sample_id)

    return guidelines


def compute_cart_ylim(
    spec: dict,
    plot_data,
    guidelines: list[dict] | dict | None = None,
) -> tuple | None:
    """Resolve y-limits while keeping every guideline visible."""

    dynamic_ylim = spec.get("cart_ylim_dynamic")

    if dynamic_ylim is None:
        cart_ylim = spec.get("cart_ylim")

    elif dynamic_ylim["mode"] == "max_plus":
        max_value = plot_data[dynamic_ylim["column"]].max()

        cart_ylim = (
            dynamic_ylim.get(
                "lower",
                0,
            ),
            max_value
            + dynamic_ylim.get(
                "offset",
                0,
            ),
        )

    else:
        raise ValueError(f"Unsupported dynamic y-limit mode: {dynamic_ylim['mode']}")

    if guidelines is None:
        return cart_ylim

    if isinstance(
        guidelines,
        dict,
    ):
        guidelines = [guidelines]

    if not guidelines:
        return cart_ylim

    y_values = plot_data.select(
        pl.col(spec["y_var"]).cast(pl.Float64, strict=False).alias("_y_value")
    )

    finite_values = (
        y_values.filter(
            pl.col("_y_value").is_not_null() & pl.col("_y_value").is_finite()
        )
        .get_column("_y_value")
        .to_list()
    )

    required_values = [
        0.0,
    ]

    required_values.extend(float(value) for value in finite_values)

    for guideline in guidelines:
        guideline_value = float(guideline["value"])

        annotation_offset = float(
            guideline.get(
                "ann_y_offset",
                0,
            )
            or 0
        )

        required_values.append(guideline_value)

        required_values.append(guideline_value + annotation_offset)

    required_lower = min(required_values)
    required_upper = max(required_values)

    if cart_ylim is None:
        lower_limit = min(
            0.0,
            required_lower,
        )
        upper_limit = required_upper
    else:
        lower_limit, upper_limit = cart_ylim

        if lower_limit is None:
            lower_limit = min(
                0.0,
                required_lower,
            )
        else:
            lower_limit = min(
                float(lower_limit),
                required_lower,
            )

        if upper_limit is None:
            upper_limit = required_upper
        else:
            upper_limit = max(
                float(upper_limit),
                required_upper,
            )

    span = max(
        upper_limit - lower_limit,
        abs(upper_limit),
        1.0,
    )

    padding = span * 0.08

    if required_upper >= upper_limit:
        upper_limit += padding

    if required_lower <= lower_limit and required_lower < 0:
        lower_limit -= padding

    return (
        lower_limit,
        upper_limit,
    )


def prepare_bar_plot_data(
    table: pl.DataFrame,
    spec: dict,
) -> pl.DataFrame:
    """Prepare one table for a bar plot."""

    filters = []

    for column_name in spec.get(
        "na_filter_columns",
        [],
    ):
        filters.append(valid_metric_expr(column_name))

    for filter_spec in spec.get(
        "filters",
        [],
    ):
        filters.append(build_filter_expression(filter_spec))

    if filters:
        combined_filter = filters[0]

        for filter_expression in filters[1:]:
            combined_filter = combined_filter & filter_expression

        table = table.filter(combined_filter)

    table = table.with_columns(
        build_value_expression(
            spec["value_spec"],
            spec["y_var"],
        )
    )

    if spec["x_var"] == "SAMPLE_ID" and {"RUN_INDEX", "SAMPLE_ID"}.issubset(
        table.columns
    ):
        max_sample_id_length = table.select(
            pl.col("SAMPLE_ID").cast(pl.Utf8).str.len_chars().max()
        ).item()

        table = table.with_columns(
            pl.when(pl.col("RUN_INDEX").is_not_null())
            .then(
                pl.col("RUN_INDEX").cast(pl.Utf8)
                + " | "
                + pl.col("SAMPLE_ID").cast(pl.Utf8).str.pad_end(max_sample_id_length)
            )
            .otherwise(pl.col("SAMPLE_ID").cast(pl.Utf8))
            .alias("PLOT_SAMPLE_ID")
        )

    if spec["fill_var"] in {"RUN", "RUN_ID"} and {
        "RUN_INDEX",
        spec["fill_var"],
    }.issubset(table.columns):
        table = table.with_columns(
            pl.when(pl.col("RUN_INDEX").is_not_null())
            .then(
                pl.col("RUN_INDEX").cast(pl.Utf8)
                + " | "
                + pl.col(spec["fill_var"]).cast(pl.Utf8)
            )
            .otherwise(pl.col(spec["fill_var"]).cast(pl.Utf8))
            .alias("PLOT_RUN")
        )

    return table


def save_plot(
    pdf_handle: PdfPages,
    plot,
) -> None:
    """Render a plotnine plot to the PDF."""

    figure = plot.draw()

    page_number = pdf_handle.get_pagecount() + 1

    figure.text(
        0.985,
        0.015,
        f"Page {page_number}",
        ha="right",
        va="bottom",
        fontsize=8,
    )

    pdf_handle.savefig(
        figure,
        bbox_inches="tight",
    )
    plt.close(figure)


def render_bar_plot(
    pdf_handle: PdfPages,
    spec: dict,
    tables: dict,
    workflow: str,
) -> None:
    """Render one bar plot."""

    plot_data = prepare_bar_plot_data(
        tables[spec["source"]],
        spec,
    )

    if spec.get("skip_if_empty") and plot_data.is_empty():
        return

    guidelines = get_available_guidelines(
        tables,
        spec,
    )

    plot_x_var = (
        "PLOT_SAMPLE_ID"
        if spec["x_var"] == "SAMPLE_ID" and "PLOT_SAMPLE_ID" in plot_data.columns
        else spec["x_var"]
    )
    plot_fill_var = (
        "PLOT_RUN"
        if spec["fill_var"] in {"RUN", "RUN_ID"} and "PLOT_RUN" in plot_data.columns
        else spec["fill_var"]
    )
    plot_x_lab = (
        "Run index | Sample ID" if plot_x_var == "PLOT_SAMPLE_ID" else spec["x_lab"]
    )

    plot = plot_bar_metric(
        data=plot_data,
        x_var=plot_x_var,
        y_var=spec["y_var"],
        fill_var=plot_fill_var,
        guide_title=spec.get(
            "guide_title",
            "Run",
        ),
        x_lab=plot_x_lab,
        y_lab=spec["y_lab"],
        cart_ylim=compute_cart_ylim(
            spec,
            plot_data,
            guidelines,
        ),
        title=resolve_plot_title(spec, workflow),
        x_lab_angle=spec.get(
            "x_lab_angle",
            ANGLE_X_NAMES,
        ),
        fig_size=spec.get(
            "fig_size",
            (15, 6),
        ),
        alpha_value=spec.get(
            "alpha_value",
            0.8,
        ),
    )
    for guideline_index, guideline in enumerate(guidelines):
        line_alpha = HV_LINE_ALPHA if guideline["alpha"] is None else guideline["alpha"]

        line_color = HV_LINE_COLOR if guideline["color"] is None else guideline["color"]

        label_x = plot_data.height + 0.9 + guideline_index * 0.6

        plot = (
            plot
            + geom_hline(
                yintercept=guideline["value"],
                alpha=line_alpha,
                color=line_color,
                size=HV_LINE_SIZE,
            )
            + annotate(
                "text",
                label_x,
                guideline["value"] + guideline["ann_y_offset"],
                label=guideline["label"],
                angle=90,
                size=12,
            )
        )

    overlay_labels = spec.get("overlay_labels")

    if overlay_labels is not None:
        plot = plot + geom_text(
            aes(label=overlay_labels["label_column"]),
            angle=overlay_labels.get(
                "angle",
                0,
            ),
            va=overlay_labels.get("va"),
            ha=overlay_labels.get("ha"),
            nudge_x=overlay_labels.get(
                "nudge_x",
                0,
            ),
            nudge_y=overlay_labels.get(
                "nudge_y",
                0,
            ),
            size=overlay_labels.get(
                "size",
                10,
            ),
        )

    save_plot(
        pdf_handle,
        plot,
    )


def render_cluster_density_scatter(
    pdf_handle: PdfPages,
    spec: dict,
    tables: dict,
    workflow: str,
) -> None:
    """Render cluster density against estimated yield."""

    plot_table = (
        tables[spec["source"]]
        .filter(valid_metric_expr("ESTIMATED_YIELD"))
        .with_columns(
            [
                pl.col("CLUSTER_DENSITY").cast(pl.Float64),
                pl.col("ESTIMATED_YIELD").cast(pl.Float64),
                (pl.col("RUN_INDEX").cast(pl.Utf8) + " | " + pl.col("RUN_ID")).alias(
                    "RUN_LABEL"
                ),
            ]
        )
    )

    if spec.get("skip_if_empty") and plot_table.is_empty():
        return

    plot = (
        ggplot(
            plot_table,
            aes(
                x="CLUSTER_DENSITY",
                y="ESTIMATED_YIELD",
            ),
        )
        + geom_point(
            aes(color="RUN_LABEL"),
            alpha=0.8,
            size=4.5,
        )
        + scale_color_manual(values=TABLEAU_20)
        + geom_text(
            aes(label="RUN_INDEX"),
            alpha=0.8,
            size=10,
            nudge_y=-5,
        )
        + guides(color=guide_legend("Run"))
        + coord_cartesian(
            xlim=(0, 400),
            ylim=(0, 150),
        )
        + xlab("Cluster density")
        + ylab("Estimated yield")
        + ggtitle(resolve_plot_title(spec, workflow))
        + scale_x_continuous(
            breaks=[index * 50 for index in range(9)],
            minor_breaks=[index * 25 for index in range(17)],
        )
        + scale_y_continuous(
            breaks=[index * 10 for index in range(16)],
            minor_breaks=[index * 5 for index in range(31)],
        )
        + theme_minimal()
        + theme(
            axis_text_x=element_text(
                angle=spec.get(
                    "x_lab_angle",
                    ANGLE_X_NAMES,
                ),
                vjust=0.5,
                hjust=0.5,
                size=7,
            ),
            figure_size=spec.get(
                "fig_size",
                (15, 6),
            ),
            axis_line=element_line(color=AXIS_LINE_COLOR),
            plot_title=element_text(
                weight="bold",
                ha="left",
            ),
            panel_grid_major=element_line(
                color=AXIS_LINE_COLOR,
                size=0.5,
                alpha=0.5,
            ),
            panel_grid_minor=element_line(
                color=AXIS_LINE_COLOR,
                size=0.25,
                alpha=0.5,
            ),
        )
    )

    save_plot(
        pdf_handle,
        plot,
    )


def render_contamination_scatter(
    pdf_handle: PdfPages,
    spec: dict,
    tables: dict,
    workflow: str,
) -> None:
    """Render a DNA contamination scatter plot."""

    plot_table = (
        tables[spec["source"]]
        .filter(
            valid_metric_expr("DNA_CONTAMINATION_SCORE")
            & valid_metric_expr("DNA_CONTAMINATION_P_VALUE")
        )
        .with_columns(
            [
                pl.col("DNA_CONTAMINATION_SCORE").cast(pl.Float64),
                pl.col("DNA_CONTAMINATION_P_VALUE").cast(pl.Float64),
                (
                    pl.col("RUN_INDEX").cast(pl.Utf8)
                    + " | "
                    + pl.col("RUN").cast(pl.Utf8)
                ).alias("RUN_LABEL"),
            ]
        )
    )

    if spec.get("skip_if_empty") and plot_table.is_empty():
        return

    score_result = (
        tables["dna_guideline_table"]
        .filter(pl.col("SAMPLE_ID") == "USL_Guideline")
        .with_columns(
            pl.when(pl.col("DNA_CONTAMINATION_SCORE") == "NA")
            .then(None)
            .otherwise(pl.col("DNA_CONTAMINATION_SCORE"))
            .cast(pl.Float64)
            .alias("DNA_CONTAMINATION_SCORE")
        )
        .select("DNA_CONTAMINATION_SCORE")
        .head(1)
        .item()
    )

    usl_contamination_score = (
        float(score_result) if score_result is not None else 5000.0
    )

    pval_result = (
        tables["dna_guideline_table"]
        .filter(pl.col("SAMPLE_ID") == "USL_Guideline")
        .with_columns(
            pl.when(pl.col("DNA_CONTAMINATION_P_VALUE") == "NA")
            .then(None)
            .otherwise(pl.col("DNA_CONTAMINATION_P_VALUE"))
            .cast(pl.Float64)
            .alias("DNA_CONTAMINATION_P_VALUE")
        )
        .select("DNA_CONTAMINATION_P_VALUE")
        .head(1)
        .item()
    )

    usl_contamination_pval = float(pval_result) if pval_result is not None else 0.05

    max_contamination_score = max(
        5000,
        plot_table.select(pl.col("DNA_CONTAMINATION_SCORE").max()).item(),
    )

    plot_color_var = "RUN_LABEL" if spec["color_var"] == "RUN" else spec["color_var"]

    plot = plot_contamination_scatter(
        data=plot_table,
        color_var=plot_color_var,
        label_var=spec["label_var"],
        guide_title=spec.get(
            "guide_title",
            "Run",
        ),
        title=resolve_plot_title(spec, workflow),
        x_lab_angle=spec.get(
            "x_lab_angle",
            ANGLE_X_NAMES,
        ),
        fig_size=spec.get(
            "fig_size",
            (15, 6),
        ),
        max_contamination_score=(max_contamination_score),
        usl_contamination_score=(usl_contamination_score),
        usl_contamination_pval=(usl_contamination_pval),
        color_values=spec.get("color_values"),
    )

    save_plot(
        pdf_handle,
        plot,
    )


def render_plot(
    pdf_handle: PdfPages,
    spec_name: str,
    spec: dict,
    tables: dict,
    workflow: str,
) -> None:
    """Dispatch one specification to its renderer."""

    plot_kind = spec["plot_kind"]

    if plot_kind == "bar":
        render_bar_plot(
            pdf_handle,
            spec,
            tables,
            workflow,
        )
        return

    if plot_kind == "cluster_density_scatter":
        render_cluster_density_scatter(pdf_handle, spec, tables, workflow)
        return

    if plot_kind == "contamination_scatter":
        render_contamination_scatter(pdf_handle, spec, tables, workflow)
        return

    raise ValueError(f"Unsupported plot kind for {spec_name}: {plot_kind}")


def build_tables(
    joint_qc_table: pl.DataFrame,
    metrics_table: pl.DataFrame,
    workflow: str,
) -> dict:
    """Build logical table subsets used by the plot specifications."""

    workflow = workflow.strip().lower()

    if workflow not in SUPPORTED_WORKFLOWS:
        message = (
            f"Unsupported workflow: {workflow}. "
            f"Expected one of: "
            f"{', '.join(sorted(SUPPORTED_WORKFLOWS))}."
        )
        logger.error(message)
        raise ValueError(message)

    metrics_table = metrics_table.sort("RUN_INDEX")
    joint_qc_table = joint_qc_table.sort("RUN_INDEX")

    metrics_table = metrics_table.filter(
        pl.col("WORKFLOW_TYPE").str.to_lowercase() == workflow
    )

    joint_qc_table = joint_qc_table.filter(
        pl.col("WORKFLOW_TYPE").str.to_lowercase() == workflow
    )

    if metrics_table.is_empty():
        message = f"No metrics rows available for {workflow} plotting."
        logger.error(message)
        raise ValueError(message)

    latest_run_index = (
        metrics_table.filter(
            pl.col("RECORD_TYPE").is_in(
                [
                    "DNA_SAMPLE",
                    "RNA_SAMPLE",
                    "SAMPLE",
                ]
            )
        )
        .select(pl.col("RUN_INDEX").cast(pl.Int64).min())
        .item()
    )

    joint_qc_guideline_table = joint_qc_table.filter(
        pl.col("RUN_ID").is_in(
            [
                "LSL_Guideline",
                "USL_Guideline",
                "Internal Guideline",
            ]
        )
    )

    guideline_table = metrics_table.filter(
        pl.col("SAMPLE_ID").is_in(
            [
                "LSL_Guideline",
                "USL_Guideline",
            ]
        )
    )

    internal_guideline_table = metrics_table.filter(
        pl.col("SAMPLE_ID") == "Internal Guideline"
    )

    dna_guideline_table = guideline_table

    rna_guideline_table = guideline_table

    data_table = metrics_table.filter(
        ~pl.col("RECORD_TYPE").is_in(
            [
                "LOWER_THRESHOLD",
                "UPPER_THRESHOLD",
            ]
        )
    )

    dna_data_table = data_table.filter(pl.col("RECORD_TYPE") == "DNA_SAMPLE")

    rna_data_table = data_table.filter(pl.col("RECORD_TYPE") == "RNA_SAMPLE")

    dna_data_table = dna_data_table.with_columns(
        [
            pl.when(pl.col("RUN_INDEX").cast(pl.Int64) == latest_run_index)
            .then(pl.lit("True"))
            .otherwise(pl.lit("False"))
            .alias("highlighted_run"),
            pl.when(pl.col("RUN_INDEX").cast(pl.Int64) == latest_run_index)
            .then(pl.col("SAMPLE_ID"))
            .otherwise(pl.lit(""))
            .alias("contamination_label"),
        ]
    )

    return {
        "joint_qc_table": joint_qc_table,
        "joint_qc_guideline_table": (joint_qc_guideline_table),
        "merged_tables": metrics_table,
        "guideline_table": guideline_table,
        "internal_guideline_table": (internal_guideline_table),
        "dna_guideline_table": (dna_guideline_table),
        "rna_guideline_table": (rna_guideline_table),
        "data_table": data_table,
        "dna_data_table": dna_data_table,
        "rna_data_table": rna_data_table,
        "dna_sample_count": (dna_data_table.height),
        "rna_sample_count": (rna_data_table.height),
    }


def generate_qc_plots(
    metrics_table: pl.DataFrame,
    joint_qc_table: pl.DataFrame,
    workflow: str,
    output_pdf: Path,
) -> None:
    """Generate a workflow-specific QC PDF."""

    workflow = workflow.strip().lower()

    if workflow not in SUPPORTED_WORKFLOWS:
        message = (
            f"Unsupported workflow: {workflow}. "
            f"Expected one of: "
            f"{', '.join(sorted(SUPPORTED_WORKFLOWS))}."
        )
        logger.error(message)
        raise ValueError(message)

    validate_plot_specs(PLOT_SPECS)

    tables = build_tables(
        joint_qc_table=joint_qc_table,
        metrics_table=metrics_table,
        workflow=workflow,
    )

    logger.info(
        f"Generating {workflow} QC plots with "
        f"{tables['dna_sample_count']} DNA sample(s) and "
        f"{tables['rna_sample_count']} RNA sample(s)."
    )

    with PdfPages(output_pdf) as pdf_handle:
        eligible_specs = [
            (spec_name, spec)
            for spec_name, spec in PLOT_SPECS.items()
            if spec.get(
                workflow,
                {},
            ).get(
                "plot",
                False,
            )
        ]

        for spec_name, spec in sorted(
            eligible_specs,
            key=lambda item: item[1][workflow]["index"],
        ):
            requires_samples = spec.get("requires_samples")

            if requires_samples == "dna" and tables["dna_sample_count"] == 0:
                continue

            if requires_samples == "rna" and tables["rna_sample_count"] == 0:
                continue

            logger.info(
                f"Rendering plot {spec_name} (index {spec[workflow]['index']})."
            )

            render_plot(
                pdf_handle,
                spec_name,
                spec,
                tables,
                workflow,
            )

    logger.info(f"QC plots written to {output_pdf}.")
