"""Unit tests for workflow-specific metric plotting helpers."""

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


def test_duplicate_plot_indices_are_rejected():
    """Duplicate positive indices within one workflow are rejected."""
    duplicate_specs = {
        "FIRST": _minimal_bar_spec(index=1),
        "SECOND": _minimal_bar_spec(index=1),
    }

    with pytest.raises(
        KeyError,
        match="duplicate plot index 1",
    ):
        _validate_plot_specs(duplicate_specs)


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


def test_missing_workflow_key_is_rejected():
    """Every plot specification must contain both workflows."""
    spec = _minimal_bar_spec()
    del spec["dragen"]

    with pytest.raises(
        KeyError,
        match="Missing workflow routing key",
    ):
        _validate_plot_specs({"TEST": spec})


def test_missing_workflow_plot_field_is_rejected():
    """Workflow routing requires the plot field."""
    spec = _minimal_bar_spec()
    spec["localapp"] = {"index": 1}

    with pytest.raises(
        KeyError,
        match="missing fields",
    ):
        _validate_plot_specs({"TEST": spec})


def test_missing_workflow_index_field_is_rejected():
    """Workflow routing requires the index field."""
    spec = _minimal_bar_spec()
    spec["localapp"] = {"plot": True}

    with pytest.raises(
        KeyError,
        match="missing fields",
    ):
        _validate_plot_specs({"TEST": spec})


def test_non_boolean_plot_flag_is_rejected():
    """Plot routing flag must be boolean."""
    spec = _minimal_bar_spec()
    spec["localapp"]["plot"] = "yes"

    with pytest.raises(
        KeyError,
        match="'plot' must be bool",
    ):
        _validate_plot_specs({"TEST": spec})


@pytest.mark.parametrize(
    "invalid_index",
    [
        -1,
        1.5,
        "1",
        None,
    ],
)
def test_invalid_plot_index_is_rejected(invalid_index):
    """Workflow plot index must be a non-negative integer."""
    spec = _minimal_bar_spec()
    spec["localapp"]["index"] = invalid_index

    with pytest.raises(
        KeyError,
        match="'index' must be",
    ):
        _validate_plot_specs({"TEST": spec})


def test_unknown_plot_kind_is_rejected():
    """Unknown plot renderer types are rejected."""
    spec = _minimal_bar_spec()
    spec["plot_kind"] = "unknown_plot"

    with pytest.raises(
        KeyError,
        match="not recognized",
    ):
        _validate_plot_specs({"TEST": spec})


def test_missing_common_plot_field_is_rejected():
    """Common required specification fields are validated."""
    spec = _minimal_bar_spec()
    del spec["source"]

    with pytest.raises(
        KeyError,
        match="Missing fields",
    ):
        _validate_plot_specs({"TEST": spec})


def test_missing_bar_specific_field_is_rejected():
    """Bar specifications require their bar-specific fields."""
    spec = _minimal_bar_spec()
    del spec["value_spec"]

    with pytest.raises(
        KeyError,
        match="value_spec",
    ):
        _validate_plot_specs({"TEST": spec})


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


def test_compute_cart_ylim_contains_both_guidelines():
    """Axis limits include all available thresholds."""

    data = pl.DataFrame(
        {
            "VALUE": [
                0.5,
                1.0,
                3.0,
            ]
        }
    )

    result = _compute_cart_ylim(
        {
            "y_var": "VALUE",
        },
        data,
        [
            {
                "value": 1.0,
                "ann_y_offset": 0,
            },
            {
                "value": 8.0,
                "ann_y_offset": 0,
            },
        ],
    )

    assert result[0] == 0
    assert result[1] > 8


# ---------------------------------------------------------------------------
# _build_filter_expression
# ---------------------------------------------------------------------------


def test_build_filter_expression_contains():
    """Contains filters perform string matching."""
    frame = pl.DataFrame(
        {
            "NAME": [
                "DNA_SAMPLE_A",
                "RNA_SAMPLE_A",
                "DNA_SAMPLE_B",
            ]
        }
    )

    result = frame.filter(
        _build_filter_expression(
            {
                "column": "NAME",
                "contains": "DNA",
            }
        )
    )

    assert result["NAME"].to_list() == [
        "DNA_SAMPLE_A",
        "DNA_SAMPLE_B",
    ]


def test_build_filter_expression_equals():
    """Equals filters keep matching values."""
    frame = pl.DataFrame(
        {
            "TYPE": [
                "DNA",
                "RNA",
                "DNA",
            ]
        }
    )

    result = frame.filter(
        _build_filter_expression(
            {
                "column": "TYPE",
                "equals": "DNA",
            }
        )
    )

    assert result["TYPE"].to_list() == [
        "DNA",
        "DNA",
    ]


def test_build_filter_expression_not_equals():
    """Not-equals filters remove matching values."""
    frame = pl.DataFrame(
        {
            "TYPE": [
                "DNA",
                "RNA",
                "SAMPLE",
            ]
        }
    )

    result = frame.filter(
        _build_filter_expression(
            {
                "column": "TYPE",
                "not_equals": "RNA",
            }
        )
    )

    assert result["TYPE"].to_list() == [
        "DNA",
        "SAMPLE",
    ]


def test_build_filter_expression_rejects_unknown_operation():
    """Unsupported filter definitions raise ValueError."""
    with pytest.raises(
        ValueError,
        match="Unsupported filter specification",
    ):
        _build_filter_expression(
            {
                "column": "TYPE",
                "startswith": "DNA",
            }
        )


# ---------------------------------------------------------------------------
# _build_value_expression
# ---------------------------------------------------------------------------


def test_build_value_expression_cast():
    """Cast operation converts values to the requested dtype."""
    frame = pl.DataFrame({"VALUE": ["1.5", "2.5"]})

    result = frame.select(
        _build_value_expression(
            {
                "operation": "cast",
                "column": "VALUE",
                "dtype": pl.Float64,
            },
            alias_name="RESULT",
        )
    )

    assert result["RESULT"].to_list() == [
        1.5,
        2.5,
    ]


def test_build_value_expression_cast_without_dtype():
    """Cast operation can return the original expression."""
    frame = pl.DataFrame({"VALUE": ["1", "2"]})

    result = frame.select(
        _build_value_expression(
            {
                "operation": "cast",
                "column": "VALUE",
            },
            alias_name="RESULT",
        )
    )

    assert result["RESULT"].to_list() == [
        "1",
        "2",
    ]


def test_build_value_expression_default_operation_is_cast():
    """Missing operation defaults to cast."""
    frame = pl.DataFrame({"VALUE": ["1", "2"]})

    result = frame.select(
        _build_value_expression(
            {
                "column": "VALUE",
                "dtype": pl.Int64,
            },
            alias_name="RESULT",
        )
    )

    assert result["RESULT"].to_list() == [
        1,
        2,
    ]


def test_build_value_expression_divide():
    """Divide operation scales the selected metric."""
    frame = pl.DataFrame({"VALUE": ["10", "20"]})

    result = frame.select(
        _build_value_expression(
            {
                "operation": "divide",
                "column": "VALUE",
                "divisor": 10,
                "dtype": pl.Float64,
            },
            alias_name="RESULT",
        )
    )

    assert result["RESULT"].to_list() == [
        1.0,
        2.0,
    ]


def test_build_value_expression_ratio():
    """Ratio operation divides numerator by denominator."""
    frame = pl.DataFrame(
        {
            "NUM": ["10", "20"],
            "DEN": ["2", "4"],
        }
    )

    result = frame.select(
        _build_value_expression(
            {
                "operation": "ratio",
                "numerator": "NUM",
                "denominator": "DEN",
                "dtype": pl.Float64,
            },
            alias_name="RESULT",
        )
    )

    assert result["RESULT"].to_list() == [
        5.0,
        5.0,
    ]


def test_build_value_expression_ratio_with_numerator_divisor():
    """Ratio supports scaling the numerator before division."""
    frame = pl.DataFrame(
        {
            "NUM": ["100", "200"],
            "DEN": ["2", "4"],
        }
    )

    result = frame.select(
        _build_value_expression(
            {
                "operation": "ratio",
                "numerator": "NUM",
                "denominator": "DEN",
                "numerator_divisor": 10,
                "dtype": pl.Float64,
            },
            alias_name="RESULT",
        )
    )

    assert result["RESULT"].to_list() == [
        5.0,
        5.0,
    ]


def test_build_value_expression_rejects_unknown_operation():
    """Unsupported value operations raise ValueError."""
    with pytest.raises(
        ValueError,
        match="Unsupported value operation",
    ):
        _build_value_expression(
            {
                "operation": "multiply",
                "column": "VALUE",
            }
        )


# ---------------------------------------------------------------------------
# _get_guideline_value
# ---------------------------------------------------------------------------


def test_get_guideline_value_returns_guideline():
    """Guideline values and plotting metadata are extracted."""
    tables = {
        "guidelines": pl.DataFrame(
            {
                "SAMPLE_ID": ["USL_Guideline"],
                "VALUE": ["12.5"],
            }
        )
    }

    result = _get_guideline_value(
        tables,
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
    )

    assert result == {
        "value": 12.5,
        "label": "USL: 12.5",
        "alpha": 0.5,
        "color": "red",
        "ann_y_offset": 2,
    }


def test_get_guideline_value_supports_explicit_label():
    """Explicit guideline labels override generated labels."""
    tables = {
        "guidelines": pl.DataFrame(
            {
                "SAMPLE_ID": ["USL_Guideline"],
                "VALUE": ["10"],
            }
        )
    }

    result = _get_guideline_value(
        tables,
        {
            "table": "guidelines",
            "sample_id": "USL_Guideline",
            "value_spec": {
                "column": "VALUE",
                "dtype": pl.Float64,
            },
            "label_prefix": "Unused",
            "label": "Custom guideline",
        },
    )

    assert result["label"] == "Custom guideline"


def test_get_guideline_value_supports_custom_id_column():
    """Guidelines can select rows using a non-SAMPLE_ID column."""
    tables = {
        "guidelines": pl.DataFrame(
            {
                "RUN_ID": ["LSL_Guideline"],
                "VALUE": ["85"],
            }
        )
    }

    result = _get_guideline_value(
        tables,
        {
            "table": "guidelines",
            "id_column": "RUN_ID",
            "sample_id": "LSL_Guideline",
            "value_spec": {
                "column": "VALUE",
                "dtype": pl.Float64,
            },
            "label_prefix": "LSL",
        },
    )

    assert result["value"] == 85.0


def test_get_guideline_value_empty_table_returns_none():
    """Empty guideline tables do not produce a guideline."""
    tables = {
        "guidelines": pl.DataFrame(
            schema={
                "SAMPLE_ID": pl.String,
                "VALUE": pl.String,
            }
        )
    }

    result = _get_guideline_value(
        tables,
        {
            "table": "guidelines",
            "sample_id": "USL_Guideline",
            "value_spec": {
                "column": "VALUE",
            },
            "label_prefix": "USL",
        },
    )

    assert result is None


def test_get_guideline_value_missing_id_column_returns_none():
    """Missing guideline ID columns are handled gracefully."""
    tables = {
        "guidelines": pl.DataFrame(
            {
                "OTHER": ["USL_Guideline"],
                "VALUE": ["10"],
            }
        )
    }

    result = _get_guideline_value(
        tables,
        {
            "table": "guidelines",
            "sample_id": "USL_Guideline",
            "value_spec": {
                "column": "VALUE",
            },
            "label_prefix": "USL",
        },
    )

    assert result is None


def test_get_guideline_value_missing_sample_returns_none():
    """Missing guideline rows return None."""
    tables = {
        "guidelines": pl.DataFrame(
            {
                "SAMPLE_ID": ["LSL_Guideline"],
                "VALUE": ["10"],
            }
        )
    }

    result = _get_guideline_value(
        tables,
        {
            "table": "guidelines",
            "sample_id": "USL_Guideline",
            "value_spec": {
                "column": "VALUE",
            },
            "label_prefix": "USL",
        },
    )

    assert result is None


def test_get_guideline_value_missing_metric_column_returns_none():
    """Missing guideline metric columns return None."""
    tables = {
        "guidelines": pl.DataFrame(
            {
                "SAMPLE_ID": ["USL_Guideline"],
            }
        )
    }

    result = _get_guideline_value(
        tables,
        {
            "table": "guidelines",
            "sample_id": "USL_Guideline",
            "value_spec": {
                "column": "MISSING",
            },
            "label_prefix": "USL",
        },
    )

    assert result is None


def test_get_guideline_value_null_metric_returns_none():
    """Null guideline metric values return None."""
    tables = {
        "guidelines": pl.DataFrame(
            {
                "SAMPLE_ID": ["USL_Guideline"],
                "VALUE": [None],
            },
            schema={
                "SAMPLE_ID": pl.String,
                "VALUE": pl.Float64,
            },
        )
    }

    result = _get_guideline_value(
        tables,
        {
            "table": "guidelines",
            "sample_id": "USL_Guideline",
            "value_spec": {
                "column": "VALUE",
            },
            "label_prefix": "USL",
        },
    )

    assert result is None


def test_get_guideline_value_literal_na_returns_none():
    """Literal NA guideline values are treated as unavailable."""

    tables = {
        "guidelines": pl.DataFrame(
            {
                "SAMPLE_ID": [
                    "USL_Guideline",
                ],
                "VALUE": [
                    "NA",
                ],
            }
        )
    }

    result = _get_guideline_value(
        tables,
        {
            "table": "guidelines",
            "sample_id": "USL_Guideline",
            "value_spec": {
                "column": "VALUE",
            },
            "python_cast": float,
            "label_prefix": "USL",
        },
    )

    assert result is None


def test_get_guideline_value_zero_lsl_returns_none():
    """A zero LSL is not treated as a drawable guideline."""

    tables = {
        "guidelines": pl.DataFrame(
            {
                "SAMPLE_ID": [
                    "LSL_Guideline",
                ],
                "VALUE": [
                    "0",
                ],
            }
        )
    }

    result = _get_guideline_value(
        tables,
        {
            "table": "guidelines",
            "sample_id": "LSL_Guideline",
            "value_spec": {
                "column": "VALUE",
            },
            "python_cast": float,
            "label_prefix": "LSL_Guideline",
        },
    )

    assert result is None


# ---------------------------------------------------------------------------
# _compute_cart_ylim
# ---------------------------------------------------------------------------


def test_compute_cart_ylim_returns_static_limits():
    """Static axis limits are returned unchanged."""
    data = pl.DataFrame({"VALUE": [1, 2, 3]})

    assert _compute_cart_ylim(
        {
            "cart_ylim": (0, 10),
        },
        data,
    ) == (0, 10)


def test_compute_cart_ylim_returns_none_without_configuration():
    """Missing axis-limit configuration returns None."""
    data = pl.DataFrame({"VALUE": [1, 2]})

    assert _compute_cart_ylim({}, data) is None


def test_compute_cart_ylim_dynamic_max_plus():
    """Dynamic limits use data maximum plus configured offset."""
    data = pl.DataFrame({"VALUE": [10, 20, 15]})

    result = _compute_cart_ylim(
        {
            "cart_ylim_dynamic": {
                "mode": "max_plus",
                "column": "VALUE",
                "offset": 5,
            }
        },
        data,
    )

    assert result == (0, 25)


def test_compute_cart_ylim_dynamic_custom_lower():
    """Dynamic limits support a custom lower bound."""
    data = pl.DataFrame({"VALUE": [10, 20]})

    result = _compute_cart_ylim(
        {
            "cart_ylim_dynamic": {
                "mode": "max_plus",
                "column": "VALUE",
                "lower": 5,
                "offset": 10,
            }
        },
        data,
    )

    assert result == (5, 30)


def test_compute_cart_ylim_rejects_unknown_dynamic_mode():
    """Unsupported dynamic limit modes raise ValueError."""
    data = pl.DataFrame({"VALUE": [1, 2]})

    with pytest.raises(
        ValueError,
        match="Unsupported dynamic y-limit mode",
    ):
        _compute_cart_ylim(
            {
                "cart_ylim_dynamic": {
                    "mode": "unknown",
                    "column": "VALUE",
                }
            },
            data,
        )


def test_compute_cart_ylim_guideline_above_bars_expands_limit():
    """Guidelines above all bars remain visible."""
    data = pl.DataFrame(
        {
            "VALUE": [
                1.0,
                2.0,
                3.0,
            ]
        }
    )

    result = _compute_cart_ylim(
        {
            "y_var": "VALUE",
        },
        data,
        {
            "value": 8.0,
            "ann_y_offset": 1.0,
        },
    )

    assert result[0] == 0
    assert result[1] > 9.0


def test_compute_cart_ylim_bar_above_guideline_uses_bar_max():
    """Bars above the guideline determine the visible upper range."""
    data = pl.DataFrame(
        {
            "VALUE": [
                3.0,
                10.0,
                15.0,
            ]
        }
    )

    result = _compute_cart_ylim(
        {
            "y_var": "VALUE",
        },
        data,
        {
            "value": 8.0,
            "ann_y_offset": 0,
        },
    )

    assert result[0] == 0
    assert result[1] > 15.0


def test_compute_cart_ylim_keeps_existing_limit_when_guideline_is_visible():
    """Configured limits remain unchanged when they already show the guideline."""
    data = pl.DataFrame(
        {
            "VALUE": [
                1.0,
                2.0,
                3.0,
            ]
        }
    )

    result = _compute_cart_ylim(
        {
            "y_var": "VALUE",
            "cart_ylim": (0, 20),
        },
        data,
        {
            "value": 8.0,
            "ann_y_offset": 1.0,
        },
    )

    assert result == (0, 20)


def test_compute_cart_ylim_expands_existing_limit_for_guideline():
    """Configured limits expand when a guideline would otherwise be clipped."""
    data = pl.DataFrame(
        {
            "VALUE": [
                1.0,
                2.0,
                3.0,
            ]
        }
    )

    result = _compute_cart_ylim(
        {
            "y_var": "VALUE",
            "cart_ylim": (0, 5),
        },
        data,
        {
            "value": 8.0,
            "ann_y_offset": 0,
        },
    )

    assert result[0] == 0
    assert result[1] > 8.0


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


def test_build_tables_rejects_unknown_workflow():
    """Unsupported workflows raise ValueError."""
    with pytest.raises(
        ValueError,
        match="Unsupported workflow",
    ):
        _build_tables(
            joint_qc_table=_joint_qc_frame(),
            metrics_table=_metrics_frame(),
            workflow="unknown",
        )


def test_build_tables_rejects_empty_selected_workflow():
    """A workflow with no matching metrics cannot be plotted."""
    with pytest.raises(
        ValueError,
        match="No metrics rows available",
    ):
        _build_tables(
            joint_qc_table=_joint_qc_frame(),
            metrics_table=_metrics_frame(),
            workflow="localapp",
        )


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
    """DRAGEN uses contamination-score bar plot at index 10."""
    spec = PLOT_SPECS["DNA_CONTAMINATION_SCORE"]

    assert spec["dragen"] == {
        "plot": True,
        "index": 10,
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


def test_render_plot_dispatches_bar(
    monkeypatch,
):
    """Bar specifications use _render_bar_plot."""
    renderer = MagicMock()

    monkeypatch.setattr(
        plotting,
        "_render_bar_plot",
        renderer,
    )

    pdf = MagicMock()

    spec = {
        "plot_kind": "bar",
    }

    _render_plot(
        pdf,
        "TEST",
        spec,
        {},
        "dragen",
    )

    renderer.assert_called_once_with(
        pdf,
        spec,
        {},
        "dragen",
    )


def test_render_plot_dispatches_cluster_density(
    monkeypatch,
):
    """Cluster-density specifications use the correct renderer."""
    renderer = MagicMock()

    monkeypatch.setattr(
        plotting,
        "_render_cluster_density_scatter",
        renderer,
    )

    pdf = MagicMock()

    spec = {
        "plot_kind": "cluster_density_scatter",
    }

    _render_plot(
        pdf,
        "TEST",
        spec,
        {},
        "dragen",
    )

    renderer.assert_called_once_with(
        pdf,
        spec,
        {},
        "dragen",
    )


def test_render_plot_dispatches_contamination(
    monkeypatch,
):
    """Contamination specifications use the correct renderer."""
    renderer = MagicMock()

    monkeypatch.setattr(
        plotting,
        "_render_contamination_scatter",
        renderer,
    )

    pdf = MagicMock()

    spec = {
        "plot_kind": "contamination_scatter",
    }

    _render_plot(
        pdf,
        "TEST",
        spec,
        {},
        "localapp",
    )

    renderer.assert_called_once_with(
        pdf,
        spec,
        {},
        "localapp",
    )


def test_render_plot_rejects_unknown_kind():
    """Unsupported renderer types raise ValueError."""
    with pytest.raises(
        ValueError,
        match="Unsupported plot kind",
    ):
        _render_plot(
            MagicMock(),
            "BROKEN",
            {
                "plot_kind": "unknown",
            },
            {},
            "dragen",
        )


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
):
    """Unknown workflow values are rejected before rendering."""
    with pytest.raises(
        ValueError,
        match="Unsupported workflow",
    ):
        Generate_qc_plots(
            metrics_table=_metrics_frame(),
            joint_qc_table=_joint_qc_frame(),
            workflow="unknown",
            output_pdf=(tmp_path / "test.pdf"),
        )


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
