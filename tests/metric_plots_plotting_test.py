"""Unit tests for workflow-specific metric plotting helpers."""

from contextlib import nullcontext
from unittest.mock import MagicMock

import polars as pl
import pytest

import tsoppy.metric_plots.plotting as plotting
from tsoppy.metric_plots.plotting import (
    Generate_qc_plots,
    _build_filter_expression,
    _build_tables,
    _build_value_expression,
    _compute_cart_ylim,
    _get_available_guidelines,
    _get_guideline_value,
    _prepare_bar_plot_data,
    _render_bar_plot,
    _render_plot,
    _resolve_plot_title,
    _save_plot,
    _valid_metric_expr,
    _validate_plot_specs,
)
from tsoppy.metric_plots.specs.plot_specs_workflows import PLOT_SPECS


def _metrics_frame(
    workflow: str = "dragen",
) -> pl.DataFrame:
    """Build synthetic metrics rows for plotting-table tests."""
    return pl.DataFrame(
        {
            "SAMPLE_ID": [
                "DNA_LATEST",
                "RNA_LATEST",
                "DNA_OLDER",
                "LSL_Guideline",
                "USL_Guideline",
                "Internal Guideline",
            ],
            "RUN": [
                "RUN_002",
                "RUN_002",
                "RUN_010",
                "GUIDELINE",
                "GUIDELINE",
                "GUIDELINE",
            ],
            "RUN_INDEX": [
                "002",
                "002",
                "010",
                "999",
                "999",
                "999",
            ],
            "WORKFLOW_TYPE": [workflow] * 6,
            "RECORD_TYPE": [
                "DNA_SAMPLE",
                "RNA_SAMPLE",
                "DNA_SAMPLE",
                "LOWER_THRESHOLD",
                "UPPER_THRESHOLD",
                "SAMPLE",
            ],
            "DNA_CONTAMINATION_SCORE": [
                "100",
                "200",
                "110",
                "0",
                "1457",
                None,
            ],
            "DNA_CONTAMINATION_P_VALUE": [
                "0.01",
                "0.02",
                "0.03",
                "0.0",
                "0.05",
                None,
            ],
            "RNA_MEDIAN_CV_GENE_500X": [
                "0.3",
                "0.2",
                "0.4",
                "0.0",
                "1.0",
                None,
            ],
            "VALUE": [
                "10",
                "20",
                "30",
                "0",
                "100",
                None,
            ],
        }
    )


def _joint_qc_frame(
    workflow: str = "dragen",
) -> pl.DataFrame:
    """Build synthetic joint-QC rows for plotting-table tests."""
    return pl.DataFrame(
        {
            "RUN_ID": [
                "RUN_002",
                "RUN_010",
                "LSL_Guideline",
                "USL_Guideline",
                "Internal Guideline",
            ],
            "RUN_INDEX": [
                "002",
                "010",
                "999",
                "999",
                "999",
            ],
            "WORKFLOW_TYPE": [workflow] * 5,
            "PCT_PF_READS": [
                "90",
                "91",
                "80",
                "100",
                "85",
            ],
            "PCT_Q30_R1": [
                "92",
                "93",
                "85",
                "100",
                "90",
            ],
            "PCT_Q30_R2": [
                "89",
                "90",
                "82",
                "100",
                "87",
            ],
            "CLUSTER_DENSITY": [
                "200",
                "210",
                "0",
                "0",
                "0",
            ],
            "ESTIMATED_YIELD": [
                "100",
                "110",
                "0",
                "0",
                "0",
            ],
            "CLUSTERS_PASSING_FILTER": [
                "95",
                "96",
                "90",
                "100",
                "92",
            ],
        }
    )


def _minimal_bar_spec(
    index: int = 1,
) -> dict:
    """Return a structurally valid minimal bar specification."""
    return {
        "localapp": {
            "plot": True,
            "index": index,
        },
        "dragen": {
            "plot": False,
            "index": 0,
        },
        "plot_kind": "bar",
        "source": "data_table",
        "title": "Test plot",
        "x_var": "SAMPLE_ID",
        "y_var": "VALUE",
        "fill_var": "RUN",
        "x_lab": "Sample ID",
        "y_lab": "Value",
        "value_spec": {
            "operation": "cast",
            "column": "VALUE",
            "dtype": pl.Float64,
        },
    }


# ---------------------------------------------------------------------------
# _resolve_plot_title
# ---------------------------------------------------------------------------


def test_resolve_plot_title_shared_title():
    """Shared plot titles are returned unchanged."""
    spec = {"title": "Shared title"}

    assert _resolve_plot_title(spec, "dragen") == "Shared title"
    assert _resolve_plot_title(spec, "localapp") == "Shared title"


def test_resolve_plot_title_workflow_specific():
    """Workflow-specific titles use the selected workflow."""
    spec = {
        "title": {
            "dragen": "Dragen title",
            "localapp": "LocalApp title",
        }
    }

    assert _resolve_plot_title(spec, "dragen") == "Dragen title"
    assert _resolve_plot_title(spec, "localapp") == "LocalApp title"


def test_workflow_specific_title_resolution_production_specs():
    """Production workflow-specific titles resolve correctly."""
    q30_spec = PLOT_SPECS["PCT_Q30_R1"]

    assert _resolve_plot_title(
        q30_spec,
        "localapp",
    ).startswith("[LocalApp")

    assert _resolve_plot_title(
        q30_spec,
        "dragen",
    ).startswith("[Dragen")


# ---------------------------------------------------------------------------
# _validate_plot_specs
# ---------------------------------------------------------------------------


def test_current_plot_specs_validate():
    """The production plot specification set is structurally valid."""
    _validate_plot_specs(PLOT_SPECS)


def test_duplicate_plot_indices_are_rejected(caplog):
    """Duplicate positive indices within one workflow are rejected."""
    duplicate_specs = {
        "FIRST": _minimal_bar_spec(index=1),
        "SECOND": _minimal_bar_spec(index=1),
    }

    with pytest.raises(KeyError):
        _validate_plot_specs(duplicate_specs)

    assert "duplicate plot index 1" in caplog.text


def test_zero_plot_indices_may_repeat():
    """Disabled plots may all use index zero."""
    first = _minimal_bar_spec(index=1)
    second = _minimal_bar_spec(index=2)

    first["dragen"] = {
        "plot": False,
        "index": 0,
    }
    second["dragen"] = {
        "plot": False,
        "index": 0,
    }

    _validate_plot_specs(
        {
            "FIRST": first,
            "SECOND": second,
        }
    )


@pytest.mark.parametrize(
    "invalid_index",
    [
        -1,
        1.5,
        "1",
        None,
    ],
)
def test_invalid_plot_index_is_rejected(invalid_index, caplog):
    """Workflow plot index must be a non-negative integer."""
    spec = _minimal_bar_spec()
    spec["localapp"]["index"] = invalid_index

    with pytest.raises(KeyError):
        _validate_plot_specs({"TEST": spec})

    assert "'index' must be" in caplog.text


def _delete_dragen_key(spec):
    del spec["dragen"]


def _set_localapp_missing_plot_field(spec):
    spec["localapp"] = {"index": 1}


def _set_localapp_missing_index_field(spec):
    spec["localapp"] = {"plot": True}


def _set_non_boolean_plot_flag(spec):
    spec["localapp"]["plot"] = "yes"


def _set_unknown_plot_kind(spec):
    spec["plot_kind"] = "unknown_plot"


def _delete_source_field(spec):
    del spec["source"]


def _delete_value_spec_field(spec):
    del spec["value_spec"]


@pytest.mark.parametrize(
    "mutate_spec, want_caplog_substring",
    [
        (
            # every plot specification must contain both workflows
            _delete_dragen_key,
            "Missing workflow routing key",
        ),
        (
            # workflow routing requires the plot field
            _set_localapp_missing_plot_field,
            "missing fields",
        ),
        (
            # workflow routing requires the index field
            _set_localapp_missing_index_field,
            "missing fields",
        ),
        (
            # plot routing flag must be boolean
            _set_non_boolean_plot_flag,
            "'plot' must be bool",
        ),
        (
            # unknown plot renderer types are rejected
            _set_unknown_plot_kind,
            "not recognized",
        ),
        (
            # common required specification fields are validated
            _delete_source_field,
            "Missing fields",
        ),
        (
            # bar specifications require their bar-specific fields
            _delete_value_spec_field,
            "value_spec",
        ),
    ],
)
def test_minimal_bar_spec_rejections(mutate_spec, want_caplog_substring, caplog):
    """_validate_plot_specs rejects a minimal bar spec that violates one rule."""
    spec = _minimal_bar_spec()
    mutate_spec(spec)

    with pytest.raises(KeyError):
        _validate_plot_specs({"TEST": spec})

    assert want_caplog_substring in caplog.text


# ---------------------------------------------------------------------------
# _valid_metric_expr
# ---------------------------------------------------------------------------


def test_valid_metric_expr_keeps_valid_values():
    """Normal metric values pass the valid-metric filter."""
    frame = pl.DataFrame(
        {
            "VALUE": [
                "1",
                "2.5",
                "0",
            ]
        }
    )

    result = frame.filter(_valid_metric_expr("VALUE"))

    assert result["VALUE"].to_list() == [
        "1",
        "2.5",
        "0",
    ]


def test_valid_metric_expr_removes_null_and_na():
    """Null and literal NA values are treated as missing."""
    frame = pl.DataFrame(
        {
            "VALUE": [
                "1",
                None,
                "NA",
                "2",
            ]
        }
    )

    result = frame.filter(_valid_metric_expr("VALUE"))

    assert result["VALUE"].to_list() == [
        "1",
        "2",
    ]


def test_get_available_guidelines_returns_lsl_and_usl():
    """Both available threshold types are returned."""

    tables = {
        "dna_guideline_table": pl.DataFrame(
            {
                "SAMPLE_ID": [
                    "LSL_Guideline",
                    "USL_Guideline",
                ],
                "VALUE": [
                    "2",
                    "8",
                ],
            }
        )
    }

    spec = {
        "source": "dna_data_table",
        "y_var": "VALUE",
        "value_spec": {
            "operation": "cast",
            "column": "VALUE",
            "dtype": pl.Float64,
        },
    }

    result = _get_available_guidelines(
        tables,
        spec,
    )

    assert len(result) == 2

    assert len({guideline["sample_id"] for guideline in result}) == len(result)

    assert [guideline["sample_id"] for guideline in result] == [
        "LSL_Guideline",
        "USL_Guideline",
    ]

    assert [guideline["value"] for guideline in result] == [
        2.0,
        8.0,
    ]


def test_get_available_guidelines_skips_na_threshold():
    """Unavailable threshold values are not plotted."""

    tables = {
        "dna_guideline_table": pl.DataFrame(
            {
                "SAMPLE_ID": [
                    "LSL_Guideline",
                    "USL_Guideline",
                ],
                "VALUE": [
                    "NA",
                    "8",
                ],
            }
        )
    }

    spec = {
        "source": "dna_data_table",
        "y_var": "VALUE",
        "value_spec": {
            "operation": "cast",
            "column": "VALUE",
            "dtype": pl.Float64,
        },
    }

    result = _get_available_guidelines(
        tables,
        spec,
    )

    assert len(result) == 1

    assert result[0]["sample_id"] == "USL_Guideline"

    assert result[0]["value"] == 8.0


# ---------------------------------------------------------------------------
# _build_filter_expression
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "inputs, exception, want",
    [
        (
            # contains: keeps rows whose column value contains the substring
            (
                "NAME",
                ["DNA_SAMPLE_A", "RNA_SAMPLE_A", "DNA_SAMPLE_B"],
                {"column": "NAME", "contains": "DNA"},
            ),
            nullcontext(),
            ["DNA_SAMPLE_A", "DNA_SAMPLE_B"],
        ),
        (
            # equals: keeps only exactly matching rows
            (
                "TYPE",
                ["DNA", "RNA", "DNA"],
                {"column": "TYPE", "equals": "DNA"},
            ),
            nullcontext(),
            ["DNA", "DNA"],
        ),
        (
            # not_equals: removes only exactly matching rows
            (
                "TYPE",
                ["DNA", "RNA", "SAMPLE"],
                {"column": "TYPE", "not_equals": "RNA"},
            ),
            nullcontext(),
            ["DNA", "SAMPLE"],
        ),
        (
            # an unsupported filter key raises ValueError
            (
                "TYPE",
                ["DNA"],
                {"column": "TYPE", "startswith": "DNA"},
            ),
            pytest.raises(ValueError),
            "Unsupported filter specification",
        ),
    ],
)
def test_build_filter_expression(inputs, exception, want, caplog):
    column, values, filter_spec = inputs
    frame = pl.DataFrame({column: values})

    with exception:
        result = frame.filter(_build_filter_expression(filter_spec))

    if isinstance(want, str):
        assert want in caplog.text
    else:
        assert result[column].to_list() == want


# ---------------------------------------------------------------------------
# _build_value_expression
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "inputs, exception, want",
    [
        (
            # cast converts values to the requested dtype
            (
                {"VALUE": ["1.5", "2.5"]},
                {"operation": "cast", "column": "VALUE", "dtype": pl.Float64},
            ),
            nullcontext(),
            [1.5, 2.5],
        ),
        (
            # cast without a dtype returns the original expression
            (
                {"VALUE": ["1", "2"]},
                {"operation": "cast", "column": "VALUE"},
            ),
            nullcontext(),
            ["1", "2"],
        ),
        (
            # a missing operation defaults to cast
            (
                {"VALUE": ["1", "2"]},
                {"column": "VALUE", "dtype": pl.Int64},
            ),
            nullcontext(),
            [1, 2],
        ),
        (
            # divide scales the selected metric
            (
                {"VALUE": ["10", "20"]},
                {
                    "operation": "divide",
                    "column": "VALUE",
                    "divisor": 10,
                    "dtype": pl.Float64,
                },
            ),
            nullcontext(),
            [1.0, 2.0],
        ),
        (
            # ratio divides numerator by denominator
            (
                {"NUM": ["10", "20"], "DEN": ["2", "4"]},
                {
                    "operation": "ratio",
                    "numerator": "NUM",
                    "denominator": "DEN",
                    "dtype": pl.Float64,
                },
            ),
            nullcontext(),
            [5.0, 5.0],
        ),
        (
            # ratio supports scaling the numerator before division
            (
                {"NUM": ["100", "200"], "DEN": ["2", "4"]},
                {
                    "operation": "ratio",
                    "numerator": "NUM",
                    "denominator": "DEN",
                    "numerator_divisor": 10,
                    "dtype": pl.Float64,
                },
            ),
            nullcontext(),
            [5.0, 5.0],
        ),
        (
            # an unsupported value operation raises ValueError
            (
                {"VALUE": ["1"]},
                {"operation": "multiply", "column": "VALUE"},
            ),
            pytest.raises(ValueError),
            "Unsupported value operation",
        ),
    ],
)
def test_build_value_expression(inputs, exception, want, caplog):
    columns, value_spec = inputs
    frame = pl.DataFrame(columns)

    with exception:
        result = frame.select(_build_value_expression(value_spec, alias_name="RESULT"))

    if isinstance(want, str):
        assert want in caplog.text
    else:
        assert result["RESULT"].to_list() == want


# ---------------------------------------------------------------------------
# _get_guideline_value
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "tables, spec, want",
    [
        (
            # guideline values and plotting metadata are extracted
            {
                "guidelines": pl.DataFrame(
                    {"SAMPLE_ID": ["USL_Guideline"], "VALUE": ["12.5"]}
                )
            },
            {
                "table": "guidelines",
                "sample_id": "USL_Guideline",
                "value_spec": {
                    "operation": "cast",
                    "column": "VALUE",
                    "dtype": pl.Float64,
                },
                "python_cast": float,
                "label_prefix": "USL",
                "alpha": 0.5,
                "color": "red",
                "ann_y_offset": 2,
            },
            {
                "value": 12.5,
                "label": "USL: 12.5",
                "alpha": 0.5,
                "color": "red",
                "ann_y_offset": 2,
            },
        ),
        (
            # explicit guideline labels override generated labels
            {
                "guidelines": pl.DataFrame(
                    {"SAMPLE_ID": ["USL_Guideline"], "VALUE": ["10"]}
                )
            },
            {
                "table": "guidelines",
                "sample_id": "USL_Guideline",
                "value_spec": {"column": "VALUE", "dtype": pl.Float64},
                "label_prefix": "Unused",
                "label": "Custom guideline",
            },
            {"label": "Custom guideline"},
        ),
        (
            # guidelines can select rows using a non-SAMPLE_ID column
            {
                "guidelines": pl.DataFrame(
                    {"RUN_ID": ["LSL_Guideline"], "VALUE": ["85"]}
                )
            },
            {
                "table": "guidelines",
                "id_column": "RUN_ID",
                "sample_id": "LSL_Guideline",
                "value_spec": {"column": "VALUE", "dtype": pl.Float64},
                "label_prefix": "LSL",
            },
            {"value": 85.0},
        ),
        (
            # empty guideline tables do not produce a guideline
            {
                "guidelines": pl.DataFrame(
                    schema={"SAMPLE_ID": pl.String, "VALUE": pl.String}
                )
            },
            {
                "table": "guidelines",
                "sample_id": "USL_Guideline",
                "value_spec": {"column": "VALUE"},
                "label_prefix": "USL",
            },
            None,
        ),
        (
            # missing guideline ID columns are handled gracefully
            {"guidelines": pl.DataFrame({"OTHER": ["USL_Guideline"], "VALUE": ["10"]})},
            {
                "table": "guidelines",
                "sample_id": "USL_Guideline",
                "value_spec": {"column": "VALUE"},
                "label_prefix": "USL",
            },
            None,
        ),
        (
            # missing guideline rows return None
            {
                "guidelines": pl.DataFrame(
                    {"SAMPLE_ID": ["LSL_Guideline"], "VALUE": ["10"]}
                )
            },
            {
                "table": "guidelines",
                "sample_id": "USL_Guideline",
                "value_spec": {"column": "VALUE"},
                "label_prefix": "USL",
            },
            None,
        ),
        (
            # missing guideline metric columns return None
            {"guidelines": pl.DataFrame({"SAMPLE_ID": ["USL_Guideline"]})},
            {
                "table": "guidelines",
                "sample_id": "USL_Guideline",
                "value_spec": {"column": "MISSING"},
                "label_prefix": "USL",
            },
            None,
        ),
        (
            # null guideline metric values return None
            {
                "guidelines": pl.DataFrame(
                    {"SAMPLE_ID": ["USL_Guideline"], "VALUE": [None]},
                    schema={"SAMPLE_ID": pl.String, "VALUE": pl.Float64},
                )
            },
            {
                "table": "guidelines",
                "sample_id": "USL_Guideline",
                "value_spec": {"column": "VALUE"},
                "label_prefix": "USL",
            },
            None,
        ),
        (
            # literal NA guideline values are treated as unavailable
            {
                "guidelines": pl.DataFrame(
                    {"SAMPLE_ID": ["USL_Guideline"], "VALUE": ["NA"]}
                )
            },
            {
                "table": "guidelines",
                "sample_id": "USL_Guideline",
                "value_spec": {"column": "VALUE"},
                "python_cast": float,
                "label_prefix": "USL",
            },
            None,
        ),
        (
            # a zero LSL is not treated as a drawable guideline
            {
                "guidelines": pl.DataFrame(
                    {"SAMPLE_ID": ["LSL_Guideline"], "VALUE": ["0"]}
                )
            },
            {
                "table": "guidelines",
                "sample_id": "LSL_Guideline",
                "value_spec": {"column": "VALUE"},
                "python_cast": float,
                "label_prefix": "LSL_Guideline",
            },
            None,
        ),
    ],
)
def test_get_guideline_value(tables, spec, want):
    result = _get_guideline_value(tables, spec)

    if want is None:
        assert result is None
    else:
        for key, value in want.items():
            assert result[key] == value


# ---------------------------------------------------------------------------
# _compute_cart_ylim
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "inputs, exception, want",
    [
        (
            # a configured static limit is returned unchanged
            ({"cart_ylim": (0, 10)}, [1, 2, 3], None),
            nullcontext(),
            (0, 10),
        ),
        (
            # missing axis-limit configuration returns None
            ({}, [1, 2], None),
            nullcontext(),
            None,
        ),
        (
            # dynamic limits use data maximum plus a configured offset
            (
                {
                    "cart_ylim_dynamic": {
                        "mode": "max_plus",
                        "column": "VALUE",
                        "offset": 5,
                    }
                },
                [10, 20, 15],
                None,
            ),
            nullcontext(),
            (0, 25),
        ),
        (
            # dynamic limits support a custom lower bound
            (
                {
                    "cart_ylim_dynamic": {
                        "mode": "max_plus",
                        "column": "VALUE",
                        "lower": 5,
                        "offset": 10,
                    }
                },
                [10, 20],
                None,
            ),
            nullcontext(),
            (5, 30),
        ),
        (
            # a configured limit is kept unchanged when the guideline already fits
            (
                {"y_var": "VALUE", "cart_ylim": (0, 20)},
                [1.0, 2.0, 3.0],
                {"value": 8.0, "ann_y_offset": 1.0},
            ),
            nullcontext(),
            (0, 20),
        ),
        (
            # an unsupported dynamic mode raises ValueError
            (
                {"cart_ylim_dynamic": {"mode": "unknown", "column": "VALUE"}},
                [1, 2],
                None,
            ),
            pytest.raises(ValueError),
            "Unsupported dynamic y-limit mode",
        ),
    ],
)
def test_compute_cart_ylim(inputs, exception, want, caplog):
    spec, data_values, guidelines = inputs
    data = pl.DataFrame({"VALUE": data_values})

    with exception:
        if guidelines is None:
            result = _compute_cart_ylim(spec, data)
        else:
            result = _compute_cart_ylim(spec, data, guidelines)

    if isinstance(want, str):
        assert want in caplog.text
    else:
        assert result == want


@pytest.mark.parametrize(
    "data_values, cart_ylim, guidelines, want_upper_exceeds",
    [
        (
            # guidelines both above and below the bars remain visible
            [0.5, 1.0, 3.0],
            None,
            [
                {"value": 1.0, "ann_y_offset": 0},
                {"value": 8.0, "ann_y_offset": 0},
            ],
            8,
        ),
        (
            # a guideline above all bars remains visible
            [1.0, 2.0, 3.0],
            None,
            {"value": 8.0, "ann_y_offset": 1.0},
            9.0,
        ),
        (
            # bars above the guideline determine the visible upper range
            [3.0, 10.0, 15.0],
            None,
            {"value": 8.0, "ann_y_offset": 0},
            15.0,
        ),
        (
            # a configured limit expands when it would otherwise clip a guideline
            [1.0, 2.0, 3.0],
            (0, 5),
            {"value": 8.0, "ann_y_offset": 0},
            8.0,
        ),
    ],
)
def test_compute_cart_ylim_expands_for_guidelines(
    data_values,
    cart_ylim,
    guidelines,
    want_upper_exceeds,
):
    """Axis limits expand to keep guidelines visible without clipping bar data."""
    data = pl.DataFrame({"VALUE": data_values})
    spec = {"y_var": "VALUE"}

    if cart_ylim is not None:
        spec["cart_ylim"] = cart_ylim

    result = _compute_cart_ylim(spec, data, guidelines)

    assert result[0] == 0
    assert result[1] > want_upper_exceeds


def test_dna_chimeric_reads_has_usl_guideline():
    """DNA chimeric-read plots use the configured USL guideline."""
    spec = PLOT_SPECS["DNA_PCT_CHIMERIC_READS"]

    guideline = spec["guideline"]

    assert guideline["table"] == "dna_guideline_table"
    assert guideline["sample_id"] == "USL_Guideline"
    assert guideline["value_spec"]["column"] == "DNA_PCT_CHIMERIC_READS"
    assert guideline["python_cast"] is float


# ---------------------------------------------------------------------------
# _prepare_bar_plot_data
# ---------------------------------------------------------------------------


def test_run_index_sample_id_generation_and_padding():
    """Sample labels include RUN_INDEX and have equal width."""
    table = pl.DataFrame(
        {
            "RUN_INDEX": ["001", "002"],
            "SAMPLE_ID": ["S1", "LONG_SAMPLE"],
            "RUN": ["RUN_A", "RUN_B"],
            "VALUE": [1.0, 2.0],
        }
    )

    plot_data = _prepare_bar_plot_data(
        table,
        _minimal_bar_spec(),
    )

    labels = plot_data["PLOT_SAMPLE_ID"].to_list()

    assert labels[0].rstrip() == "001 | S1"
    assert labels[1].rstrip() == "002 | LONG_SAMPLE"
    assert len(labels[0]) == len(labels[1])


def test_prepare_bar_plot_data_preserves_original_sample_id():
    """Plot labels do not modify canonical SAMPLE_ID."""
    table = pl.DataFrame(
        {
            "RUN_INDEX": ["001"],
            "SAMPLE_ID": ["SAMPLE_A"],
            "RUN": ["RUN_A"],
            "VALUE": ["10"],
        }
    )

    result = _prepare_bar_plot_data(
        table,
        _minimal_bar_spec(),
    )

    assert result["SAMPLE_ID"].to_list() == ["SAMPLE_A"]
    assert result["PLOT_SAMPLE_ID"].to_list() == [
        "001 | SAMPLE_A",
    ]


def test_prepare_bar_plot_data_creates_run_legend():
    """Run legend labels include RUN_INDEX and RUN."""
    table = pl.DataFrame(
        {
            "RUN_INDEX": ["001", "002"],
            "SAMPLE_ID": ["S1", "S2"],
            "RUN": ["RUN_A", "RUN_B"],
            "VALUE": [1.0, 2.0],
        }
    )

    result = _prepare_bar_plot_data(
        table,
        _minimal_bar_spec(),
    )

    assert result["PLOT_RUN"].to_list() == [
        "001 | RUN_A",
        "002 | RUN_B",
    ]


def test_run_index_run_id_legend_generation():
    """Run legend labels include RUN_INDEX before RUN_ID."""
    table = pl.DataFrame(
        {
            "RUN_INDEX": ["001", "002"],
            "RUN_ID": ["RUN_A", "RUN_B"],
            "VALUE": [1.0, 2.0],
        }
    )

    spec = {
        "x_var": "RUN_ID",
        "y_var": "VALUE",
        "fill_var": "RUN_ID",
        "value_spec": {
            "operation": "cast",
            "column": "VALUE",
            "dtype": pl.Float64,
        },
    }

    plot_data = _prepare_bar_plot_data(
        table,
        spec,
    )

    assert plot_data["PLOT_RUN"].to_list() == [
        "001 | RUN_A",
        "002 | RUN_B",
    ]


def test_prepare_bar_plot_data_casts_value_column():
    """Configured value transformation is applied."""
    table = pl.DataFrame(
        {
            "SAMPLE_ID": ["S1", "S2"],
            "RUN": ["R1", "R1"],
            "VALUE": ["1.5", "2.5"],
        }
    )

    result = _prepare_bar_plot_data(
        table,
        _minimal_bar_spec(),
    )

    assert result["VALUE"].to_list() == [
        1.5,
        2.5,
    ]


def test_prepare_bar_plot_data_filters_na_values():
    """NA-filter columns remove null and literal NA rows."""
    spec = _minimal_bar_spec()
    spec["na_filter_columns"] = ["VALUE"]

    table = pl.DataFrame(
        {
            "SAMPLE_ID": ["S1", "S2", "S3", "S4"],
            "RUN": ["R1"] * 4,
            "VALUE": ["1", "NA", None, "4"],
        }
    )

    result = _prepare_bar_plot_data(
        table,
        spec,
    )

    assert result["SAMPLE_ID"].to_list() == [
        "S1",
        "S4",
    ]


def test_prepare_bar_plot_data_applies_multiple_filters():
    """Multiple configured filters are combined with AND."""
    spec = _minimal_bar_spec()
    spec["filters"] = [
        {
            "column": "TYPE",
            "equals": "DNA",
        },
        {
            "column": "SAMPLE_ID",
            "contains": "KEEP",
        },
    ]

    table = pl.DataFrame(
        {
            "SAMPLE_ID": [
                "KEEP_A",
                "DROP_A",
                "KEEP_B",
            ],
            "TYPE": [
                "DNA",
                "DNA",
                "RNA",
            ],
            "RUN": ["R1"] * 3,
            "VALUE": ["1", "2", "3"],
        }
    )

    result = _prepare_bar_plot_data(
        table,
        spec,
    )

    assert result["SAMPLE_ID"].to_list() == [
        "KEEP_A",
    ]


def test_prepare_bar_plot_data_handles_null_run_index():
    """Rows without a run index retain their original display labels."""
    table = pl.DataFrame(
        {
            "RUN_INDEX": [None, "001"],
            "SAMPLE_ID": ["GUIDELINE", "SAMPLE"],
            "RUN": ["GUIDELINE", "RUN_A"],
            "VALUE": [1.0, 2.0],
        },
        schema_overrides={
            "RUN_INDEX": pl.String,
        },
    )

    result = _prepare_bar_plot_data(
        table,
        _minimal_bar_spec(),
    )

    assert result.get_column("PLOT_SAMPLE_ID")[0].rstrip() == "GUIDELINE"

    assert result.get_column("PLOT_SAMPLE_ID")[1].rstrip() == "001 | SAMPLE"

    assert result.get_column("PLOT_RUN")[0] == "GUIDELINE"

    assert result.get_column("PLOT_RUN")[1] == "001 | RUN_A"


# ---------------------------------------------------------------------------
# _build_tables
# ---------------------------------------------------------------------------


def test_latest_run_uses_minimum_numeric_run_index():
    """Latest highlighted run is minimum numeric RUN_INDEX."""
    tables = _build_tables(
        joint_qc_table=_joint_qc_frame(),
        metrics_table=_metrics_frame(),
        workflow="dragen",
    )

    dna_table = tables["dna_data_table"]

    assert (
        dna_table.filter(pl.col("SAMPLE_ID") == "DNA_LATEST")
        .select("highlighted_run")
        .item()
        == "True"
    )

    assert (
        dna_table.filter(pl.col("SAMPLE_ID") == "DNA_OLDER")
        .select("highlighted_run")
        .item()
        == "False"
    )


def test_latest_run_comparison_is_numeric_not_lexical():
    """RUN_INDEX comparison uses numeric ordering."""
    metrics = _metrics_frame().with_columns(
        pl.when(pl.col("SAMPLE_ID") == "DNA_LATEST")
        .then(pl.lit("10"))
        .when(pl.col("SAMPLE_ID") == "DNA_OLDER")
        .then(pl.lit("2"))
        .otherwise(pl.col("RUN_INDEX"))
        .alias("RUN_INDEX")
    )

    tables = _build_tables(
        joint_qc_table=_joint_qc_frame(),
        metrics_table=metrics,
        workflow="dragen",
    )

    dna_table = tables["dna_data_table"]

    assert (
        dna_table.filter(pl.col("SAMPLE_ID") == "DNA_OLDER")
        .select("highlighted_run")
        .item()
        == "True"
    )


def test_record_type_controls_dna_rna_selection():
    """DNA/RNA plotting tables are selected from RECORD_TYPE."""
    tables = _build_tables(
        joint_qc_table=_joint_qc_frame(),
        metrics_table=_metrics_frame(),
        workflow="dragen",
    )

    assert set(tables["dna_data_table"]["SAMPLE_ID"].to_list()) == {
        "DNA_LATEST",
        "DNA_OLDER",
    }

    assert set(tables["rna_data_table"]["SAMPLE_ID"].to_list()) == {
        "RNA_LATEST",
    }


def test_build_tables_filters_workflow_case_insensitively():
    """Workflow matching is case-insensitive."""
    metrics = pl.concat(
        [
            _metrics_frame("dragen"),
            _metrics_frame("localapp"),
        ]
    )

    joint = pl.concat(
        [
            _joint_qc_frame("dragen"),
            _joint_qc_frame("localapp"),
        ]
    )

    tables = _build_tables(
        joint_qc_table=joint,
        metrics_table=metrics,
        workflow="DRAGEN",
    )

    assert set(
        tables["merged_tables"]["WORKFLOW_TYPE"].str.to_lowercase().to_list()
    ) == {"dragen"}

    assert set(
        tables["joint_qc_table"]["WORKFLOW_TYPE"].str.to_lowercase().to_list()
    ) == {"dragen"}


def test_build_tables_strips_workflow_whitespace():
    """Workflow input is normalized before matching."""
    tables = _build_tables(
        joint_qc_table=_joint_qc_frame(),
        metrics_table=_metrics_frame(),
        workflow="  DRAGEN  ",
    )

    assert tables["dna_sample_count"] == 2


def test_build_tables_rejects_unknown_workflow(caplog):
    """Unsupported workflows raise ValueError."""
    with pytest.raises(ValueError):
        _build_tables(
            joint_qc_table=_joint_qc_frame(),
            metrics_table=_metrics_frame(),
            workflow="unknown",
        )

    assert "Unsupported workflow" in caplog.text


def test_build_tables_rejects_empty_selected_workflow(caplog):
    """A workflow with no matching metrics cannot be plotted."""
    with pytest.raises(ValueError):
        _build_tables(
            joint_qc_table=_joint_qc_frame(),
            metrics_table=_metrics_frame(),
            workflow="localapp",
        )

    assert "No metrics rows available" in caplog.text


def test_build_tables_extracts_threshold_guidelines():
    """LSL and USL threshold rows are separated."""
    tables = _build_tables(
        joint_qc_table=_joint_qc_frame(),
        metrics_table=_metrics_frame(),
        workflow="dragen",
    )

    assert set(tables["guideline_table"]["SAMPLE_ID"].to_list()) == {
        "LSL_Guideline",
        "USL_Guideline",
    }


def test_build_tables_extracts_internal_guideline():
    """Internal guideline row is stored separately."""
    tables = _build_tables(
        joint_qc_table=_joint_qc_frame(),
        metrics_table=_metrics_frame(),
        workflow="dragen",
    )

    assert tables["internal_guideline_table"]["SAMPLE_ID"].to_list() == [
        "Internal Guideline"
    ]


def test_build_tables_extracts_joint_qc_guidelines():
    """Run-level guideline rows are separated."""
    tables = _build_tables(
        joint_qc_table=_joint_qc_frame(),
        metrics_table=_metrics_frame(),
        workflow="dragen",
    )

    assert set(tables["joint_qc_guideline_table"]["RUN_ID"].to_list()) == {
        "LSL_Guideline",
        "USL_Guideline",
        "Internal Guideline",
    }


def test_build_tables_removes_threshold_rows_from_data_table():
    """Threshold rows are excluded from regular sample data."""
    tables = _build_tables(
        joint_qc_table=_joint_qc_frame(),
        metrics_table=_metrics_frame(),
        workflow="dragen",
    )

    sample_ids = set(tables["data_table"]["SAMPLE_ID"].to_list())

    assert "LSL_Guideline" not in sample_ids
    assert "USL_Guideline" not in sample_ids


def test_build_tables_sets_contamination_label_only_for_latest_run():
    """Only latest-run DNA samples receive contamination labels."""
    tables = _build_tables(
        joint_qc_table=_joint_qc_frame(),
        metrics_table=_metrics_frame(),
        workflow="dragen",
    )

    dna = tables["dna_data_table"]

    assert (
        dna.filter(pl.col("SAMPLE_ID") == "DNA_LATEST")
        .select("contamination_label")
        .item()
        == "DNA_LATEST"
    )

    assert (
        dna.filter(pl.col("SAMPLE_ID") == "DNA_OLDER")
        .select("contamination_label")
        .item()
        == ""
    )


def test_build_tables_returns_sample_counts():
    """DNA and RNA sample counts reflect selected rows."""
    tables = _build_tables(
        joint_qc_table=_joint_qc_frame(),
        metrics_table=_metrics_frame(),
        workflow="dragen",
    )

    assert tables["dna_sample_count"] == 2
    assert tables["rna_sample_count"] == 1


def test_build_tables_sorts_metrics_by_run_index():
    """Metrics are sorted by RUN_INDEX."""
    tables = _build_tables(
        joint_qc_table=_joint_qc_frame(),
        metrics_table=_metrics_frame(),
        workflow="dragen",
    )

    indexes = tables["merged_tables"]["RUN_INDEX"].to_list()

    assert indexes == sorted(indexes)


# ---------------------------------------------------------------------------
# Production workflow specifications
# ---------------------------------------------------------------------------


def test_dragen_contamination_scatter_is_disabled():
    """LocalApp contamination scatters remain disabled for DRAGEN."""
    for spec_name in (
        "DNA_CONTAMINATION_P_VALUE_BY_RUN",
        "DNA_CONTAMINATION_P_VALUE_HIGHLIGHTED_RUN",
    ):
        assert PLOT_SPECS[spec_name]["dragen"] == {
            "plot": False,
            "index": 0,
        }


def test_dragen_contamination_score_bar_is_enabled():
    """DRAGEN uses contamination-score bar plot at index 20."""
    spec = PLOT_SPECS["DNA_CONTAMINATION_SCORE"]

    assert spec["dragen"] == {
        "plot": True,
        "index": 20,
    }

    assert spec["localapp"] == {
        "plot": False,
        "index": 0,
    }

    assert spec["plot_kind"] == "bar"

    assert spec["y_var"] == "DNA_CONTAMINATION_SCORE"


def test_dragen_run_level_plots_are_enabled_in_order():
    """Five shared run-level plots remain DRAGEN indices 1-5."""
    expected = {
        "CLUSTERS_PASSING_FILTER": 1,
        "ESTIMATED_YIELD": 2,
        "PCT_PF_READS": 3,
        "PCT_Q30_R1": 4,
        "PCT_Q30_R2": 5,
    }

    for (
        spec_name,
        expected_index,
    ) in expected.items():
        assert PLOT_SPECS[spec_name]["dragen"] == {
            "plot": True,
            "index": expected_index,
        }


def test_all_enabled_plot_indices_are_unique():
    """Enabled plot indices are unique within each workflow."""
    for workflow in (
        "dragen",
        "localapp",
    ):
        indices = [
            spec[workflow]["index"]
            for spec in PLOT_SPECS.values()
            if spec[workflow]["plot"]
        ]

        assert len(indices) == len(set(indices))


def test_all_enabled_plot_indices_are_positive():
    """Enabled plots always have positive indices."""
    for workflow in (
        "dragen",
        "localapp",
    ):
        for spec in PLOT_SPECS.values():
            if spec[workflow]["plot"]:
                assert spec[workflow]["index"] > 0


def test_all_disabled_plot_indices_are_zero():
    """Disabled plots consistently use index zero."""
    for workflow in (
        "dragen",
        "localapp",
    ):
        for spec in PLOT_SPECS.values():
            if not spec[workflow]["plot"]:
                assert spec[workflow]["index"] == 0


# ---------------------------------------------------------------------------
# _save_plot
# ---------------------------------------------------------------------------


def test_save_plot_draws_saves_numbers_and_closes(
    monkeypatch,
):
    """_save_plot draws, numbers, saves and closes the figure."""
    figure = MagicMock()

    plot = MagicMock()
    plot.draw.return_value = figure

    pdf_handle = MagicMock()
    pdf_handle.get_pagecount.return_value = 2

    close_mock = MagicMock()

    monkeypatch.setattr(
        plotting.plt,
        "close",
        close_mock,
    )

    _save_plot(
        pdf_handle,
        plot,
    )

    plot.draw.assert_called_once_with()

    figure.text.assert_called_once_with(
        0.985,
        0.015,
        "Page 3",
        ha="right",
        va="bottom",
        fontsize=8,
    )

    pdf_handle.savefig.assert_called_once_with(
        figure,
        bbox_inches="tight",
    )

    close_mock.assert_called_once_with(figure)


def test_save_plot_first_page_number_is_one(
    monkeypatch,
):
    """An empty PDF begins numbering at page one."""
    figure = MagicMock()

    plot = MagicMock()
    plot.draw.return_value = figure

    pdf_handle = MagicMock()
    pdf_handle.get_pagecount.return_value = 0

    monkeypatch.setattr(
        plotting.plt,
        "close",
        MagicMock(),
    )

    _save_plot(
        pdf_handle,
        plot,
    )

    assert figure.text.call_args.args[2] == "Page 1"


# ---------------------------------------------------------------------------
# _render_plot dispatch
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "plot_kind, renderer_name, workflow",
    [
        ("bar", "_render_bar_plot", "dragen"),
        ("cluster_density_scatter", "_render_cluster_density_scatter", "dragen"),
        ("contamination_scatter", "_render_contamination_scatter", "localapp"),
    ],
)
def test_render_plot_dispatches(plot_kind, renderer_name, workflow, monkeypatch):
    """Each supported plot_kind dispatches to its corresponding renderer."""
    renderer = MagicMock()
    monkeypatch.setattr(plotting, renderer_name, renderer)

    pdf = MagicMock()
    spec = {"plot_kind": plot_kind}

    _render_plot(pdf, "TEST", spec, {}, workflow)

    renderer.assert_called_once_with(pdf, spec, {}, workflow)


def test_render_plot_rejects_unknown_kind(caplog):
    """Unsupported renderer types raise ValueError."""
    with pytest.raises(ValueError):
        _render_plot(
            MagicMock(),
            "BROKEN",
            {
                "plot_kind": "unknown",
            },
            {},
            "dragen",
        )

    assert "Unsupported plot kind" in caplog.text


# ---------------------------------------------------------------------------
# _render_bar_plot
# ---------------------------------------------------------------------------


def test_render_bar_plot_calls_plot_function_and_save(
    monkeypatch,
):
    """Bar rendering passes prepared data to the plotting helper."""
    frame = pl.DataFrame(
        {
            "SAMPLE_ID": ["S1"],
            "PLOT_SAMPLE_ID": ["001 | S1"],
            "RUN": ["RUN_A"],
            "PLOT_RUN": ["001 | RUN_A"],
            "VALUE": [10.0],
        }
    )

    monkeypatch.setattr(
        plotting,
        "_prepare_bar_plot_data",
        MagicMock(return_value=frame),
    )

    plot_object = MagicMock()

    bar_mock = MagicMock(return_value=plot_object)

    monkeypatch.setattr(
        plotting,
        "Plot_bar_metric",
        bar_mock,
    )

    save_mock = MagicMock()

    monkeypatch.setattr(
        plotting,
        "_save_plot",
        save_mock,
    )

    spec = _minimal_bar_spec()

    _render_bar_plot(
        MagicMock(),
        spec,
        {
            "data_table": pl.DataFrame(),
        },
        "localapp",
    )

    kwargs = bar_mock.call_args.kwargs

    assert kwargs["x_var"] == ("PLOT_SAMPLE_ID")

    assert kwargs["fill_var"] == ("PLOT_RUN")

    assert kwargs["x_lab"] == ("Run index | Sample ID")

    assert kwargs["title"] == ("Test plot")

    save_mock.assert_called_once()


def test_render_bar_plot_skip_if_empty(
    monkeypatch,
):
    """Empty prepared bar data is skipped when configured."""
    monkeypatch.setattr(
        plotting,
        "_prepare_bar_plot_data",
        MagicMock(return_value=pl.DataFrame()),
    )

    bar_mock = MagicMock()

    monkeypatch.setattr(
        plotting,
        "Plot_bar_metric",
        bar_mock,
    )

    save_mock = MagicMock()

    monkeypatch.setattr(
        plotting,
        "_save_plot",
        save_mock,
    )

    spec = _minimal_bar_spec()
    spec["skip_if_empty"] = True

    _render_bar_plot(
        MagicMock(),
        spec,
        {
            "data_table": pl.DataFrame(),
        },
        "localapp",
    )

    bar_mock.assert_not_called()
    save_mock.assert_not_called()


def test_render_bar_plot_draws_all_available_guidelines(
    monkeypatch,
):
    """Every available LSL and USL is drawn on the bar plot."""

    frame = pl.DataFrame(
        {
            "SAMPLE_ID": ["S1"],
            "VALUE": [3.0],
            "RUN": ["RUN_A"],
        }
    )

    monkeypatch.setattr(
        plotting,
        "_prepare_bar_plot_data",
        MagicMock(return_value=frame),
    )

    guidelines = [
        {
            "sample_id": "LSL_Guideline",
            "value": 1.0,
            "label": "LSL_Guideline: 1.0",
            "alpha": 0.3,
            "color": "red",
            "ann_y_offset": 0,
        },
        {
            "sample_id": "USL_Guideline",
            "value": 8.0,
            "label": "USL_Guideline: 8.0",
            "alpha": 0.3,
            "color": "red",
            "ann_y_offset": 0,
        },
    ]

    monkeypatch.setattr(
        plotting,
        "_get_available_guidelines",
        MagicMock(return_value=guidelines),
    )

    bar_mock = MagicMock(return_value=MagicMock())

    monkeypatch.setattr(
        plotting,
        "Plot_bar_metric",
        bar_mock,
    )

    hline_mock = MagicMock(
        side_effect=[
            MagicMock(),
            MagicMock(),
        ]
    )

    monkeypatch.setattr(
        plotting,
        "geom_hline",
        hline_mock,
    )

    annotate_mock = MagicMock(
        side_effect=[
            MagicMock(),
            MagicMock(),
        ]
    )

    monkeypatch.setattr(
        plotting,
        "annotate",
        annotate_mock,
    )

    save_mock = MagicMock()

    monkeypatch.setattr(
        plotting,
        "_save_plot",
        save_mock,
    )

    spec = _minimal_bar_spec()
    label_x_positions = [call.args[1] for call in annotate_mock.call_args_list]

    assert len(set(label_x_positions)) == len(label_x_positions)
    _render_bar_plot(
        MagicMock(),
        spec,
        {
            "data_table": pl.DataFrame(),
        },
        "localapp",
    )

    assert hline_mock.call_count == 2

    line_values = [call.kwargs["yintercept"] for call in hline_mock.call_args_list]

    assert line_values == [
        1.0,
        8.0,
    ]

    labels = [call.kwargs["label"] for call in annotate_mock.call_args_list]

    assert labels == [
        "LSL_Guideline: 1.0",
        "USL_Guideline: 8.0",
    ]

    kwargs = bar_mock.call_args.kwargs

    assert kwargs["cart_ylim"][0] == 0
    assert kwargs["cart_ylim"][1] > 8.0

    save_mock.assert_called_once()


# ---------------------------------------------------------------------------
# Generate_qc_plots orchestration
# ---------------------------------------------------------------------------


def test_generate_qc_plots_rejects_unknown_workflow(
    tmp_path,
    caplog,
):
    """Unknown workflow values are rejected before rendering."""
    with pytest.raises(ValueError):
        Generate_qc_plots(
            metrics_table=_metrics_frame(),
            joint_qc_table=_joint_qc_frame(),
            workflow="unknown",
            output_pdf=(tmp_path / "test.pdf"),
        )

    assert "Unsupported workflow" in caplog.text


def test_generate_qc_plots_normalizes_workflow(
    monkeypatch,
    tmp_path,
):
    """Workflow names are stripped and lower-cased."""
    validate_mock = MagicMock()

    monkeypatch.setattr(
        plotting,
        "_validate_plot_specs",
        validate_mock,
    )

    build_mock = MagicMock(
        return_value={
            "dna_sample_count": 0,
            "rna_sample_count": 0,
        }
    )

    monkeypatch.setattr(
        plotting,
        "_build_tables",
        build_mock,
    )

    monkeypatch.setattr(
        plotting,
        "PLOT_SPECS",
        {},
    )

    fake_pdf = MagicMock()
    fake_pdf.__enter__.return_value = MagicMock()

    monkeypatch.setattr(
        plotting,
        "PdfPages",
        MagicMock(return_value=fake_pdf),
    )

    Generate_qc_plots(
        metrics_table=_metrics_frame(),
        joint_qc_table=_joint_qc_frame(),
        workflow="  DRAGEN  ",
        output_pdf=(tmp_path / "test.pdf"),
    )

    assert build_mock.call_args.kwargs["workflow"] == "dragen"


def test_generate_qc_plots_renders_enabled_specs_in_index_order(
    monkeypatch,
    tmp_path,
):
    """Enabled plots are rendered in configured index order."""
    monkeypatch.setattr(
        plotting,
        "_validate_plot_specs",
        MagicMock(),
    )

    monkeypatch.setattr(
        plotting,
        "_build_tables",
        MagicMock(
            return_value={
                "dna_sample_count": 1,
                "rna_sample_count": 1,
            }
        ),
    )

    test_specs = {
        "THIRD": {
            "dragen": {
                "plot": True,
                "index": 3,
            },
        },
        "FIRST": {
            "dragen": {
                "plot": True,
                "index": 1,
            },
        },
        "SECOND": {
            "dragen": {
                "plot": True,
                "index": 2,
            },
        },
        "DISABLED": {
            "dragen": {
                "plot": False,
                "index": 0,
            },
        },
    }

    monkeypatch.setattr(
        plotting,
        "PLOT_SPECS",
        test_specs,
    )

    rendered = []

    def fake_render(
        pdf_handle,
        spec_name,
        spec,
        tables,
        workflow,
    ):
        rendered.append(spec_name)

    monkeypatch.setattr(
        plotting,
        "_render_plot",
        fake_render,
    )

    fake_pdf = MagicMock()
    fake_pdf.__enter__.return_value = MagicMock()

    monkeypatch.setattr(
        plotting,
        "PdfPages",
        MagicMock(return_value=fake_pdf),
    )

    Generate_qc_plots(
        metrics_table=_metrics_frame(),
        joint_qc_table=_joint_qc_frame(),
        workflow="dragen",
        output_pdf=(tmp_path / "test.pdf"),
    )

    assert rendered == [
        "FIRST",
        "SECOND",
        "THIRD",
    ]


def test_generate_qc_plots_skips_dna_plot_without_dna_samples(
    monkeypatch,
    tmp_path,
):
    """DNA-only plots are skipped when no DNA samples exist."""
    monkeypatch.setattr(
        plotting,
        "_validate_plot_specs",
        MagicMock(),
    )

    monkeypatch.setattr(
        plotting,
        "_build_tables",
        MagicMock(
            return_value={
                "dna_sample_count": 0,
                "rna_sample_count": 1,
            }
        ),
    )

    monkeypatch.setattr(
        plotting,
        "PLOT_SPECS",
        {
            "DNA": {
                "dragen": {
                    "plot": True,
                    "index": 1,
                },
                "requires_samples": "dna",
            },
            "RNA": {
                "dragen": {
                    "plot": True,
                    "index": 2,
                },
                "requires_samples": "rna",
            },
        },
    )

    render_mock = MagicMock()

    monkeypatch.setattr(
        plotting,
        "_render_plot",
        render_mock,
    )

    fake_pdf = MagicMock()
    fake_pdf.__enter__.return_value = MagicMock()

    monkeypatch.setattr(
        plotting,
        "PdfPages",
        MagicMock(return_value=fake_pdf),
    )

    Generate_qc_plots(
        metrics_table=_metrics_frame(),
        joint_qc_table=_joint_qc_frame(),
        workflow="dragen",
        output_pdf=(tmp_path / "test.pdf"),
    )

    assert render_mock.call_count == 1

    assert render_mock.call_args.args[1] == "RNA"


def test_generate_qc_plots_skips_rna_plot_without_rna_samples(
    monkeypatch,
    tmp_path,
):
    """RNA-only plots are skipped when no RNA samples exist."""
    monkeypatch.setattr(
        plotting,
        "_validate_plot_specs",
        MagicMock(),
    )

    monkeypatch.setattr(
        plotting,
        "_build_tables",
        MagicMock(
            return_value={
                "dna_sample_count": 1,
                "rna_sample_count": 0,
            }
        ),
    )

    monkeypatch.setattr(
        plotting,
        "PLOT_SPECS",
        {
            "DNA": {
                "dragen": {
                    "plot": True,
                    "index": 1,
                },
                "requires_samples": "dna",
            },
            "RNA": {
                "dragen": {
                    "plot": True,
                    "index": 2,
                },
                "requires_samples": "rna",
            },
        },
    )

    render_mock = MagicMock()

    monkeypatch.setattr(
        plotting,
        "_render_plot",
        render_mock,
    )

    fake_pdf = MagicMock()
    fake_pdf.__enter__.return_value = MagicMock()

    monkeypatch.setattr(
        plotting,
        "PdfPages",
        MagicMock(return_value=fake_pdf),
    )

    Generate_qc_plots(
        metrics_table=_metrics_frame(),
        joint_qc_table=_joint_qc_frame(),
        workflow="dragen",
        output_pdf=(tmp_path / "test.pdf"),
    )

    assert render_mock.call_count == 1

    assert render_mock.call_args.args[1] == "DNA"


def test_generate_qc_plots_renders_run_plot_without_samples(
    monkeypatch,
    tmp_path,
):
    """Run-level plots do not require DNA or RNA samples."""
    monkeypatch.setattr(
        plotting,
        "_validate_plot_specs",
        MagicMock(),
    )

    monkeypatch.setattr(
        plotting,
        "_build_tables",
        MagicMock(
            return_value={
                "dna_sample_count": 0,
                "rna_sample_count": 0,
            }
        ),
    )

    monkeypatch.setattr(
        plotting,
        "PLOT_SPECS",
        {
            "RUN_METRIC": {
                "dragen": {
                    "plot": True,
                    "index": 1,
                },
                "requires_samples": None,
            },
        },
    )

    render_mock = MagicMock()

    monkeypatch.setattr(
        plotting,
        "_render_plot",
        render_mock,
    )

    fake_pdf = MagicMock()
    fake_pdf.__enter__.return_value = MagicMock()

    monkeypatch.setattr(
        plotting,
        "PdfPages",
        MagicMock(return_value=fake_pdf),
    )

    Generate_qc_plots(
        metrics_table=_metrics_frame(),
        joint_qc_table=_joint_qc_frame(),
        workflow="dragen",
        output_pdf=(tmp_path / "test.pdf"),
    )

    render_mock.assert_called_once()
