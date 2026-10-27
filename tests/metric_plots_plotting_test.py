"""Unit tests for workflow-specific metric plotting helpers."""

from contextlib import nullcontext
from unittest.mock import MagicMock

import polars as pl
from pytest import mark, raises
from matplotlib.backends.backend_pdf import PdfPages

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
    _render_cluster_density_scatter,
    _render_contamination_scatter,
    _render_plot,
    _resolve_plot_title,
    _valid_metric_expr,
    _validate_plot_specs,
)
from tsoppy.metric_plots.specs.main import PLOT_SPECS

# Test data representing bar plot spec
test_bar_plot_spec_base = {
    "title": "Test plot",
    "x_var": "SAMPLE_ID",
    "y_var": "VALUE",
    "fill_var": "RUN",
    "x_lab": "Sample ID",
    "y_lab": "Value",
}

test_bar_plot_spec_localapp = {
    "plot": True,
    "index": 1,
}

test_bar_plot_spec_dragen = {
    "plot": False,
    "index": 0,
}

test_bar_plot_spec_value_spec = {
    "operation": "cast",
    "column": "VALUE",
    "dtype": pl.Float64,
}

test_bar_plot_spec = {
    "localapp": test_bar_plot_spec_localapp,
    "dragen": test_bar_plot_spec_dragen,
    "plot_kind": "bar",
    "source": "joint_qc_table",
    "requires_samples": None,
    "title": "End-to-end QC test",
    "x_var": "RUN_ID",
    "y_var": "VALUE",
    "fill_var": "RUN_ID",
    "x_lab": "Run ID",
    "y_lab": "Percentage",
    "na_filter_columns": [
        "VALUE",
    ],
    "cart_ylim": (
        0,
        100,
    ),
    "skip_if_empty": True,
    "value_spec": test_bar_plot_spec_value_spec,
    **test_bar_plot_spec_base,
}

test_bar_plot_spec_complete = {
    "TEST_RUN_METRIC": {
        "localapp": test_bar_plot_spec_localapp,
        "dragen": test_bar_plot_spec_localapp,
        "plot_kind": "bar",
        "source": "joint_qc_table",
        "requires_samples": None,
        "title": "End-to-end QC test",
        "x_var": "RUN_ID",
        "y_var": "PCT_PF_READS",
        "fill_var": "RUN_ID",
        "x_lab": "Run ID",
        "y_lab": "Percentage",
        "na_filter_columns": [
            "PCT_PF_READS",
        ],
        "cart_ylim": (
            0,
            100,
        ),
        "skip_if_empty": True,
        "value_spec": {
            "operation": "cast",
            "column": "PCT_PF_READS",
            "dtype": pl.Float64,
        },
    }
}

# test data representing joint qc and metrics frame
test_run_id_column_guideline = [
    "LSL_Guideline",
    "USL_Guideline",
    "Internal Guideline",
]

test_run_column_guideline = ["GUIDELINE"] * 3

test_run_index_column_guideline = ["999"] * 3

test_joint_qc_frame_pct_pf_reads_column_guideline = [
    "80",
    "100",
    "85",
]

test_joint_qc_frame_pct_q30_r1_column_guideline = [
    "85",
    "100",
    "90",
]

test_joint_qc_frame_pct_q30_r2_column_guideline = [
    "82",
    "100",
    "87",
]

test_joint_qc_frame_clusters_passing_filter_column_guideline = [
    "90",
    "100",
    "92",
]

test_metrics_frame_RECORD_TYPE_column_guideline = [
    "LOWER_THRESHOLD",
    "UPPER_THRESHOLD",
    "SAMPLE",
]

test_metrics_frame_DNA_CONTAMINATION_SCORE_column_guideline = [
    "0",
    "1457",
    None,
]

test_metrics_frame_DNA_CONTAMINATION_P_VALUE_column_guideline = [
    "0.0",
    "0.05",
    None,
]

test_metrics_frame_RNA_MEDIAN_CV_GENE_500X_column_guideline = [
    "0.0",
    "1.0",
    None,
]

test_metrics_frame_VALUE_column_guideline = [
    "0",
    "100",
    None,
]

test_joint_qc_frame_guideline_table = pl.DataFrame(
    {
        "RUN_ID": test_run_id_column_guideline,
        "RUN_INDEX": test_run_index_column_guideline,
        "WORKFLOW_TYPE": ["dragen"] * 3,
        "PCT_PF_READS": test_joint_qc_frame_pct_pf_reads_column_guideline,
        "PCT_Q30_R1": test_joint_qc_frame_pct_q30_r1_column_guideline,
        "PCT_Q30_R2": test_joint_qc_frame_pct_q30_r2_column_guideline,
        "CLUSTER_DENSITY": ["0"] * 3,
        "ESTIMATED_YIELD": ["0"] * 3,
        "CLUSTERS_PASSING_FILTER": test_joint_qc_frame_clusters_passing_filter_column_guideline,
    }
)

test_metrics_frame_guideline_table = pl.DataFrame(
    {
        "SAMPLE_ID": test_run_id_column_guideline,
        "RUN": test_run_column_guideline,
        "RUN_INDEX": test_run_index_column_guideline,
        "WORKFLOW_TYPE": ["dragen"] * 3,
        "RECORD_TYPE": test_metrics_frame_RECORD_TYPE_column_guideline,
        "DNA_CONTAMINATION_SCORE": test_metrics_frame_DNA_CONTAMINATION_SCORE_column_guideline,
        "DNA_CONTAMINATION_P_VALUE": test_metrics_frame_DNA_CONTAMINATION_P_VALUE_column_guideline,
        "RNA_MEDIAN_CV_GENE_500X": test_metrics_frame_RNA_MEDIAN_CV_GENE_500X_column_guideline,
        "VALUE": test_metrics_frame_VALUE_column_guideline,
    }
)

test_joint_qc_frame_table = pl.DataFrame(
    {
        "RUN_ID": ["RUN_002", "RUN_010"] + test_run_id_column_guideline,
        "RUN_INDEX": ["002", "010"] + test_run_index_column_guideline,
        "WORKFLOW_TYPE": ["dragen"] * 5,
        "PCT_PF_READS": ["90", "91"]
        + test_joint_qc_frame_pct_pf_reads_column_guideline,
        "PCT_Q30_R1": ["92", "93"] + test_joint_qc_frame_pct_q30_r1_column_guideline,
        "PCT_Q30_R2": ["89", "90"] + test_joint_qc_frame_pct_q30_r2_column_guideline,
        "CLUSTER_DENSITY": ["200", "210"] + ["0"] * 3,
        "ESTIMATED_YIELD": ["110", "100"] + ["0"] * 3,
        "CLUSTERS_PASSING_FILTER": [
            "95",
            "96",
        ]
        + test_joint_qc_frame_clusters_passing_filter_column_guideline,
    }
)

test_metrics_frame_table = pl.DataFrame(
    {
        "SAMPLE_ID": [
            "DNA_LATEST",
            "RNA_LATEST",
            "DNA_OLDER",
        ]
        + test_run_id_column_guideline,
        "RUN": [
            "RUN_002",
            "RUN_002",
            "RUN_010",
        ]
        + test_run_column_guideline,
        "RUN_INDEX": [
            "002",
            "002",
            "010",
        ]
        + test_run_index_column_guideline,
        "WORKFLOW_TYPE": ["dragen"] * 6,
        "RECORD_TYPE": [
            "DNA_SAMPLE",
            "RNA_SAMPLE",
            "DNA_SAMPLE",
        ]
        + test_metrics_frame_RECORD_TYPE_column_guideline,
        "DNA_CONTAMINATION_SCORE": [
            "100",
            "200",
            "110",
        ]
        + test_metrics_frame_DNA_CONTAMINATION_SCORE_column_guideline,
        "DNA_CONTAMINATION_P_VALUE": [
            "0.01",
            "0.02",
            "0.03",
        ]
        + test_metrics_frame_DNA_CONTAMINATION_P_VALUE_column_guideline,
        "RNA_MEDIAN_CV_GENE_500X": [
            "0.3",
            "0.2",
            "0.4",
        ]
        + test_metrics_frame_RNA_MEDIAN_CV_GENE_500X_column_guideline,
        "VALUE": [
            "10",
            "20",
            "30",
        ]
        + test_metrics_frame_VALUE_column_guideline,
    }
)


@mark.parametrize(
    "spec, workflow, want",
    [
        (
            # a plain string title is shared across both workflows
            {"title": "Shared title"},
            "dragen",
            "Shared title",
        ),
        (
            {"title": "Shared title"},
            "localapp",
            "Shared title",
        ),
        (
            # a per-workflow title dict picks the title for the selected workflow
            {"title": {"dragen": "Dragen title", "localapp": "LocalApp title"}},
            "dragen",
            "Dragen title",
        ),
        (
            {"title": {"dragen": "Dragen title", "localapp": "LocalApp title"}},
            "localapp",
            "LocalApp title",
        ),
    ],
)
def test_resolve_plot_title(spec, workflow, want):
    assert _resolve_plot_title(spec, workflow) == want


@mark.parametrize(
    "input, exception, want",
    [
        (
            PLOT_SPECS,
            nullcontext(),
            (
                None,
                "",
            ),
        ),
        (
            {
                "FIRST": test_bar_plot_spec,
                "SECOND": test_bar_plot_spec,
            },
            raises(KeyError),
            (
                None,
                "duplicate plot index 1",
            ),
        ),
        (
            {
                "FIRST": test_bar_plot_spec,
                "SECOND": {
                    "localapp": {
                        "plot": True,
                        "index": 2,
                    },
                    "dragen": test_bar_plot_spec_dragen,
                    "plot_kind": "bar",
                    "source": "data_table",
                    "value_spec": test_bar_plot_spec_value_spec,
                    **test_bar_plot_spec_base,
                },
            },
            nullcontext(),
            (
                None,
                "",
            ),
        ),
        (
            {
                "TEST": {
                    "localapp": {
                        "plot": True,
                        "index": -1,
                    },
                    "dragen": test_bar_plot_spec_dragen,
                    "plot_kind": "bar",
                    "source": "data_table",
                    "value_spec": test_bar_plot_spec_value_spec,
                    **test_bar_plot_spec_base,
                }
            },
            raises(KeyError),
            (
                None,
                "'index' must be",
            ),
        ),
        (
            {
                "TEST": {
                    "localapp": {
                        "plot": True,
                        "index": 1.5,
                    },
                    "dragen": test_bar_plot_spec_dragen,
                    "plot_kind": "bar",
                    "source": "data_table",
                    "value_spec": test_bar_plot_spec_value_spec,
                    **test_bar_plot_spec_base,
                }
            },
            raises(KeyError),
            (
                None,
                "'index' must be",
            ),
        ),
        (
            {
                "TEST": {
                    "localapp": {
                        "plot": True,
                        "index": "1",
                    },
                    "dragen": test_bar_plot_spec_dragen,
                    "plot_kind": "bar",
                    "source": "data_table",
                    "value_spec": test_bar_plot_spec_value_spec,
                    **test_bar_plot_spec_base,
                }
            },
            raises(KeyError),
            (
                None,
                "'index' must be",
            ),
        ),
        (
            {
                "TEST": {
                    "localapp": {
                        "plot": True,
                        "index": None,
                    },
                    "dragen": test_bar_plot_spec_dragen,
                    "plot_kind": "bar",
                    "source": "data_table",
                    "value_spec": test_bar_plot_spec_value_spec,
                    **test_bar_plot_spec_base,
                }
            },
            raises(KeyError),
            (
                None,
                "'index' must be",
            ),
        ),
        (
            {
                "TEST": {
                    "localapp": {
                        "index": 1,
                    },
                    "dragen": test_bar_plot_spec_dragen,
                    "plot_kind": "bar",
                    "source": "data_table",
                    "value_spec": test_bar_plot_spec_value_spec,
                    **test_bar_plot_spec_base,
                }
            },
            raises(KeyError),
            (
                None,
                "missing fields",
            ),
        ),
        (
            {
                "TEST": {
                    "localapp": {
                        "plot": True,
                    },
                    "dragen": test_bar_plot_spec_dragen,
                    "plot_kind": "bar",
                    "source": "data_table",
                    "value_spec": test_bar_plot_spec_value_spec,
                    **test_bar_plot_spec_base,
                }
            },
            raises(KeyError),
            (
                None,
                "missing fields",
            ),
        ),
        (
            {
                "TEST": {
                    "localapp": {
                        "index": 1,
                        "plot": "yes",
                    },
                    "dragen": test_bar_plot_spec_dragen,
                    "plot_kind": "bar",
                    "source": "data_table",
                    "value_spec": test_bar_plot_spec_value_spec,
                    **test_bar_plot_spec_base,
                }
            },
            raises(KeyError),
            (
                None,
                "'plot' must be bool",
            ),
        ),
        (
            {
                "TEST": {
                    "localapp": test_bar_plot_spec_localapp,
                    "plot_kind": "bar",
                    "source": "data_table",
                    "value_spec": test_bar_plot_spec_value_spec,
                    **test_bar_plot_spec_base,
                }
            },
            raises(KeyError),
            (
                None,
                "Missing workflow routing key",
            ),
        ),
        (
            {
                "TEST": {
                    "localapp": test_bar_plot_spec_localapp,
                    "plot_kind": "bar",
                    "value_spec": test_bar_plot_spec_value_spec,
                    **test_bar_plot_spec_base,
                }
            },
            raises(KeyError),
            (
                None,
                "Missing fields",
            ),
        ),
        (
            {
                "TEST": {
                    "localapp": test_bar_plot_spec_localapp,
                    "dragen": test_bar_plot_spec_dragen,
                    "plot_kind": "bar",
                    "source": "data_table",
                    **test_bar_plot_spec_base,
                }
            },
            raises(KeyError),
            (
                None,
                "value_spec",
            ),
        ),
        (
            {
                "TEST": {
                    "localapp": test_bar_plot_spec_localapp,
                    "dragen": test_bar_plot_spec_dragen,
                    "plot_kind": "unknown_plot",
                    "source": "data_table",
                    "value_spec": test_bar_plot_spec_value_spec,
                    **test_bar_plot_spec_base,
                }
            },
            raises(KeyError),
            (
                None,
                "not recognized",
            ),
        ),
    ],
)
def test_validate_plot_specs(input, exception, want, caplog):
    with exception:
        got = _validate_plot_specs(input)
        assert got == want[0]
        assert want[1] in caplog.text


@mark.parametrize(
    "input, want",
    [
        (
            # normal metric values pass the valid-metric filter
            ["1", "2.5", "0"],
            ["1", "2.5", "0"],
        ),
        (
            # null and literal NA values are treated as missing
            ["1", None, "NA", "2"],
            ["1", "2"],
        ),
    ],
)
def test_valid_metric_expr(input, want):
    got = pl.DataFrame({"VALUE": input}).filter(_valid_metric_expr("VALUE"))
    assert got["VALUE"].to_list() == want


@mark.parametrize(
    "input, wants",
    [
        (
            # both available threshold types are returned
            ["2", "8"],
            (
                [
                    "LSL_Guideline",
                    "USL_Guideline",
                ],
                [2.0, 8.0],
            ),
        ),
        (
            # unavailable threshold values are not plotted
            ["NA", "8"],
            (
                ["USL_Guideline"],
                [8.0],
            ),
        ),
    ],
)
def test_get_available_guidelines(input, wants):
    tables = {
        "dna_guideline_table": pl.DataFrame(
            {
                "SAMPLE_ID": ["LSL_Guideline", "USL_Guideline"],
                "VALUE": input,
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

    got = _get_available_guidelines(tables, spec)

    assert len(got) == len(wants[0])
    assert len({guideline["sample_id"] for guideline in got}) == len(wants[0])
    assert [guideline["sample_id"] for guideline in got] == wants[0]
    assert [guideline["value"] for guideline in got] == wants[1]


@mark.parametrize(
    "inputs, exception, wants",
    [
        (
            # contains: keeps rows whose column value contains the substring
            (
                pl.DataFrame(
                    {"NAME": ["DNA_SAMPLE_A", "RNA_SAMPLE_A", "DNA_SAMPLE_B"]}
                ),
                {"column": "NAME", "contains": "DNA"},
            ),
            nullcontext(),
            (["DNA_SAMPLE_A", "DNA_SAMPLE_B"], ""),
        ),
        (
            # equals: keeps only exactly matching rows
            (
                pl.DataFrame({"TYPE": ["DNA", "RNA", "DNA"]}),
                {"column": "TYPE", "equals": "DNA"},
            ),
            nullcontext(),
            (["DNA", "DNA"], ""),
        ),
        (
            # equals: keeps only exactly matching rows
            (
                pl.DataFrame({"TYPE": ["DNA", "RNA", "DNA"]}),
                {"column": "TYPE", "equals": "DNA"},
            ),
            nullcontext(),
            (["DNA", "DNA"], ""),
        ),
        (
            # not_equals: removes only exactly matching rows
            (
                pl.DataFrame({"TYPE": ["DNA", "RNA", "SAMPLE"]}),
                {"column": "TYPE", "not_equals": "RNA"},
            ),
            nullcontext(),
            (["DNA", "SAMPLE"], ""),
        ),
        (
            # an unsupported filter key raises ValueError
            (
                pl.DataFrame({"TYPE": ["DNA"]}),
                {"column": "TYPE", "startswith": "DNA"},
            ),
            raises(ValueError),
            (
                None,
                "Unsupported filter specification",
            ),
        ),
    ],
)
def test_build_filter_expression(inputs, exception, wants, caplog):
    with exception:
        got = inputs[0].filter(_build_filter_expression(inputs[1]))
        column = inputs[0].columns[0]
        assert got[column].to_list() == wants[0]

    assert wants[1] in caplog.text


@mark.parametrize(
    "inputs, exception, wants",
    [
        (
            # cast converts values to the requested dtype
            (
                pl.DataFrame({"VALUE": ["1.5", "2.5"]}),
                {"operation": "cast", "column": "VALUE", "dtype": pl.Float64},
            ),
            nullcontext(),
            ([1.5, 2.5], ""),
        ),
        (
            # cast without a dtype returns the original expression
            (
                pl.DataFrame({"VALUE": ["1", "2"]}),
                {"operation": "cast", "column": "VALUE"},
            ),
            nullcontext(),
            (["1", "2"], ""),
        ),
        (
            # a missing operation defaults to cast
            (
                pl.DataFrame({"VALUE": ["1", "2"]}),
                {"column": "VALUE", "dtype": pl.Int64},
            ),
            nullcontext(),
            ([1, 2], ""),
        ),
        (
            # divide scales the selected metric
            (
                pl.DataFrame({"VALUE": ["10", "20"]}),
                {
                    "operation": "divide",
                    "column": "VALUE",
                    "divisor": 10,
                    "dtype": pl.Float64,
                },
            ),
            nullcontext(),
            ([1.0, 2.0], ""),
        ),
        (
            # ratio divides numerator by denominator
            (
                pl.DataFrame({"NUM": ["10", "20"], "DEN": ["2", "4"]}),
                {
                    "operation": "ratio",
                    "numerator": "NUM",
                    "denominator": "DEN",
                    "dtype": pl.Float64,
                },
            ),
            nullcontext(),
            ([5.0, 5.0], ""),
        ),
        (
            # ratio supports scaling the numerator before division
            (
                pl.DataFrame({"NUM": ["100", "200"], "DEN": ["2", "4"]}),
                {
                    "operation": "ratio",
                    "numerator": "NUM",
                    "denominator": "DEN",
                    "numerator_divisor": 10,
                    "dtype": pl.Float64,
                },
            ),
            nullcontext(),
            ([5.0, 5.0], ""),
        ),
        (
            # an unsupported value operation raises ValueError
            (
                pl.DataFrame({"VALUE": ["1"]}),
                {"operation": "multiply", "column": "VALUE"},
            ),
            raises(ValueError),
            (
                None,
                "Unsupported value operation",
            ),
        ),
    ],
)
def test_build_value_expression(inputs, exception, wants, caplog):
    with exception:
        got = inputs[0].select(_build_value_expression(inputs[1], alias_name="RESULT"))
        assert got["RESULT"].to_list() == wants[0]

    assert wants[1] in caplog.text


@mark.parametrize(
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
    got = _get_guideline_value(tables, spec)

    if want is None:
        assert got is want
    else:
        for key, value in want.items():
            assert got[key] == value


@mark.parametrize(
    "inputs, exception, wants",
    [
        (
            # a configured static limit is returned unchanged
            ({"cart_ylim": (0, 10)}, pl.DataFrame({"VALUE": [1, 2, 3]}), None),
            nullcontext(),
            ((0, 10), ""),
        ),
        (
            # missing axis-limit configuration returns None
            ({}, pl.DataFrame({"VALUE": [1, 2]}), None),
            nullcontext(),
            (None, ""),
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
                pl.DataFrame({"VALUE": [10, 20, 15]}),
                None,
            ),
            nullcontext(),
            ((0, 25), ""),
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
                pl.DataFrame({"VALUE": [10, 20]}),
                None,
            ),
            nullcontext(),
            ((5, 30), ""),
        ),
        (
            # a configured limit is kept unchanged when the guideline already fits
            (
                {"y_var": "VALUE", "cart_ylim": (0, 20)},
                pl.DataFrame({"VALUE": [1.0, 2.0, 3.0]}),
                {"value": 8.0, "ann_y_offset": 1.0},
            ),
            nullcontext(),
            ((0, 20), ""),
        ),
        (
            # an unsupported dynamic mode raises ValueError
            (
                {"cart_ylim_dynamic": {"mode": "unknown", "column": "VALUE"}},
                pl.DataFrame({"VALUE": [1, 2]}),
                None,
            ),
            raises(ValueError),
            (None, "Unsupported dynamic y-limit mode"),
        ),
        (
            # guidelines both above and below the bars remain visible
            (
                {"y_var": "VALUE"},
                pl.DataFrame({"VALUE": [0.5, 1.0, 3.0]}),
                [
                    {"value": 1.0, "ann_y_offset": 0},
                    {"value": 8.0, "ann_y_offset": 0},
                ],
            ),
            nullcontext(),
            ((0, 8.64), ""),
        ),
        (
            # a guideline above all bars remains visible
            (
                {"y_var": "VALUE"},
                pl.DataFrame({"VALUE": [1.0, 2.0, 3.0]}),
                {"value": 8.0, "ann_y_offset": 1.0},
            ),
            nullcontext(),
            ((0, 9.72), ""),
        ),
        (
            # bars above the guideline determine the visible upper range
            (
                {"y_var": "VALUE"},
                pl.DataFrame({"VALUE": [3.0, 10.0, 15.0]}),
                {"value": 8.0, "ann_y_offset": 0},
            ),
            nullcontext(),
            ((0, 16.2), ""),
        ),
        (
            # a configured limit expands when it would otherwise clip a guideline
            (
                {"y_var": "VALUE", "cart_ylim": (0, 5)},
                pl.DataFrame({"VALUE": [1.0, 2.0, 3.0]}),
                {"value": 8.0, "ann_y_offset": 0},
            ),
            nullcontext(),
            ((0, 8.64), ""),
        ),
    ],
)
def test_compute_cart_ylim(inputs, exception, wants, caplog):
    with exception:
        got = _compute_cart_ylim(inputs[0], inputs[1], inputs[2])
        assert got == wants[0]

    assert wants[1] in caplog.text


@mark.parametrize(
    "table, spec, want",
    [
        (
            pl.DataFrame(
                {
                    "RUN_INDEX": ["001", "002"],
                    "SAMPLE_ID": ["S1", "LONG_SAMPLE"],
                    "RUN": ["RUN_A", "RUN_B"],
                    "VALUE": [1.0, 2.0],
                }
            ),
            test_bar_plot_spec,
            {
                "expected_sample_ids": ["001 | S1         ", "002 | LONG_SAMPLE"],
                "expected_run_labels": ["001 | RUN_A", "002 | RUN_B"],
                "expected_value_rows": ["S1", "LONG_SAMPLE"],
            },
        ),
        (
            pl.DataFrame(
                {
                    "RUN_INDEX": ["001", "002"],
                    "SAMPLE_ID": ["S1", "S2"],
                    "RUN_ID": ["RUN_A", "RUN_B"],
                    "VALUE": [1.0, 2.0],
                }
            ),
            {**test_bar_plot_spec, "fill_var": "RUN_ID"},
            {
                "expected_sample_ids": ["001 | S1", "002 | S2"],
                "expected_run_labels": ["001 | RUN_A", "002 | RUN_B"],
                "expected_value_rows": ["S1", "S2"],
            },
        ),
        (
            pl.DataFrame(
                {
                    "RUN_INDEX": [None, "001"],
                    "SAMPLE_ID": ["GUIDELINE", "SAMPLE"],
                    "RUN": ["GUIDELINE", "RUN_A"],
                    "VALUE": [1.0, 2.0],
                },
                schema_overrides={
                    "RUN_INDEX": pl.String,
                },
            ),
            test_bar_plot_spec,
            {
                "expected_sample_ids": ["GUIDELINE", "001 | SAMPLE   "],
                "expected_run_labels": ["GUIDELINE", "001 | RUN_A"],
                "expected_value_rows": ["GUIDELINE", "SAMPLE"],
            },
        ),
        (
            pl.DataFrame(
                {
                    "SAMPLE_ID": ["S1", "S2", "S3", "S4"],
                    "RUN": ["R1"] * 4,
                    "VALUE": ["1", "NA", None, "4"],
                }
            ),
            {**test_bar_plot_spec, "na_filter_columns": ["VALUE"]},
            {
                "expected_sample_ids": ["S1", "S4"],
                "expected_run_labels": ["R1", "R1"],
                "expected_value_rows": ["S1", "S4"],
            },
        ),
        (
            pl.DataFrame(
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
            ),
            {
                **test_bar_plot_spec,
                "filters": [
                    {
                        "column": "TYPE",
                        "equals": "DNA",
                    },
                    {
                        "column": "SAMPLE_ID",
                        "contains": "KEEP",
                    },
                ],
            },
            {
                "expected_sample_ids": ["KEEP_A"],
                "expected_run_labels": ["R1"],
                "expected_value_rows": ["KEEP_A"],
            },
        ),
    ],
)
def test_prepare_bar_plot_data(table, spec, want):
    """Prepare bar plot data adds display labels and applies the configured filters."""
    result = _prepare_bar_plot_data(table, spec)

    if "PLOT_SAMPLE_ID" in result.columns:
        print(result["PLOT_SAMPLE_ID"].to_list())
        assert result["PLOT_SAMPLE_ID"].to_list() == want["expected_sample_ids"]

    if "PLOT_RUN" in result.columns:
        assert result["PLOT_RUN"].to_list() == want["expected_run_labels"]

    assert result["SAMPLE_ID"].to_list() == want["expected_value_rows"]


@mark.parametrize(
    "inputs, exception, wants",
    [
        # Latest highlighted run is minimum numeric RUN_INDEX.
        (
            (test_joint_qc_frame_table, test_metrics_frame_table, "dragen"),
            nullcontext(),
            (
                {
                    "joint_qc_table": test_joint_qc_frame_table,
                    "joint_qc_guideline_table": test_joint_qc_frame_guideline_table,
                    "merged_tables": test_metrics_frame_table,
                    "guideline_table": test_metrics_frame_guideline_table.slice(0, 2),
                    "internal_guideline_table": test_metrics_frame_guideline_table.slice(
                        2, 1
                    ),
                    "dna_guideline_table": test_metrics_frame_guideline_table.slice(
                        0, 2
                    ),
                    "rna_guideline_table": test_metrics_frame_guideline_table.slice(
                        0, 2
                    ),
                    "data_table": pl.concat(
                        [
                            test_metrics_frame_table.slice(0, 3),
                            test_metrics_frame_guideline_table.slice(2, 1),
                        ]
                    ),
                    "dna_data_table": pl.concat(
                        [
                            test_metrics_frame_table.slice(0, 1),
                            test_metrics_frame_table.slice(2, 1),
                        ]
                    )
                    .insert_column(9, pl.Series("highlighted_run", ["True", "False"]))
                    .insert_column(
                        10, pl.Series("contamination_label", ["DNA_LATEST", ""])
                    ),
                    "rna_data_table": test_metrics_frame_table.slice(1, 1),
                    "dna_sample_count": 2,
                    "rna_sample_count": 1,
                },
                "",
            ),
        ),
        # RUN_INDEX comparison uses numeric ordering.
        (
            (
                test_joint_qc_frame_table,
                test_metrics_frame_table.clone().replace_column(
                    2,
                    pl.Series(
                        "RUN_INDEX",
                        ["010", "010", "002"] + test_run_index_column_guideline,
                    ),
                ),
                "dragen",
            ),
            nullcontext(),
            (
                {
                    "joint_qc_table": test_joint_qc_frame_table,
                    "joint_qc_guideline_table": test_joint_qc_frame_guideline_table,
                    "merged_tables": test_metrics_frame_table.clone()
                    .replace_column(
                        2,
                        pl.Series(
                            "RUN_INDEX",
                            ["010", "010", "002"] + test_run_index_column_guideline,
                        ),
                    )
                    .sort("RUN_INDEX"),
                    "guideline_table": test_metrics_frame_guideline_table.slice(0, 2),
                    "internal_guideline_table": test_metrics_frame_guideline_table.slice(
                        2, 1
                    ),
                    "dna_guideline_table": test_metrics_frame_guideline_table.slice(
                        0, 2
                    ),
                    "rna_guideline_table": test_metrics_frame_guideline_table.slice(
                        0, 2
                    ),
                    "data_table": pl.concat(
                        [
                            test_metrics_frame_table.slice(0, 3).replace_column(
                                2, pl.Series("RUN_INDEX", ["010", "010", "002"])
                            ),
                            test_metrics_frame_guideline_table.slice(2, 1),
                        ]
                    ).sort("RUN_INDEX"),
                    "dna_data_table": pl.concat(
                        [
                            test_metrics_frame_table.slice(0, 1),
                            test_metrics_frame_table.slice(2, 1),
                        ]
                    )
                    .replace_column(2, pl.Series("RUN_INDEX", ["010", "002"]))
                    .sort("RUN_INDEX")
                    .insert_column(9, pl.Series("highlighted_run", ["True", "False"]))
                    .insert_column(
                        10, pl.Series("contamination_label", ["DNA_OLDER", ""])
                    ),
                    "rna_data_table": test_metrics_frame_table.slice(1, 1)
                    .replace_column(2, pl.Series("RUN_INDEX", ["010"]))
                    .sort("RUN_INDEX"),
                    "dna_sample_count": 2,
                    "rna_sample_count": 1,
                },
                "",
            ),
        ),
        (
            (test_joint_qc_frame_table, test_metrics_frame_table, "DRAGEN"),
            nullcontext(),
            (
                {
                    "joint_qc_table": test_joint_qc_frame_table,
                    "joint_qc_guideline_table": test_joint_qc_frame_guideline_table,
                    "merged_tables": test_metrics_frame_table,
                    "guideline_table": test_metrics_frame_guideline_table.slice(0, 2),
                    "internal_guideline_table": test_metrics_frame_guideline_table.slice(
                        2, 1
                    ),
                    "dna_guideline_table": test_metrics_frame_guideline_table.slice(
                        0, 2
                    ),
                    "rna_guideline_table": test_metrics_frame_guideline_table.slice(
                        0, 2
                    ),
                    "data_table": pl.concat(
                        [
                            test_metrics_frame_table.slice(0, 3),
                            test_metrics_frame_guideline_table.slice(2, 1),
                        ]
                    ),
                    "dna_data_table": pl.concat(
                        [
                            test_metrics_frame_table.slice(0, 1),
                            test_metrics_frame_table.slice(2, 1),
                        ]
                    )
                    .insert_column(9, pl.Series("highlighted_run", ["True", "False"]))
                    .insert_column(
                        10, pl.Series("contamination_label", ["DNA_LATEST", ""])
                    ),
                    "rna_data_table": test_metrics_frame_table.slice(1, 1),
                    "dna_sample_count": 2,
                    "rna_sample_count": 1,
                },
                "",
            ),
        ),
        (
            (test_joint_qc_frame_table, test_metrics_frame_table, "  dragen  "),
            nullcontext(),
            (
                {
                    "joint_qc_table": test_joint_qc_frame_table,
                    "joint_qc_guideline_table": test_joint_qc_frame_guideline_table,
                    "merged_tables": test_metrics_frame_table,
                    "guideline_table": test_metrics_frame_guideline_table.slice(0, 2),
                    "internal_guideline_table": test_metrics_frame_guideline_table.slice(
                        2, 1
                    ),
                    "dna_guideline_table": test_metrics_frame_guideline_table.slice(
                        0, 2
                    ),
                    "rna_guideline_table": test_metrics_frame_guideline_table.slice(
                        0, 2
                    ),
                    "data_table": pl.concat(
                        [
                            test_metrics_frame_table.slice(0, 3),
                            test_metrics_frame_guideline_table.slice(2, 1),
                        ]
                    ),
                    "dna_data_table": pl.concat(
                        [
                            test_metrics_frame_table.slice(0, 1),
                            test_metrics_frame_table.slice(2, 1),
                        ]
                    )
                    .insert_column(9, pl.Series("highlighted_run", ["True", "False"]))
                    .insert_column(
                        10, pl.Series("contamination_label", ["DNA_LATEST", ""])
                    ),
                    "rna_data_table": test_metrics_frame_table.slice(1, 1),
                    "dna_sample_count": 2,
                    "rna_sample_count": 1,
                },
                "",
            ),
        ),
        (
            (test_joint_qc_frame_table, test_metrics_frame_table, "unknown"),
            raises(ValueError),
            (
                None,
                "Unsupported workflow",
            ),
        ),
        (
            (test_joint_qc_frame_table, test_metrics_frame_table, "localapp"),
            raises(ValueError),
            (
                None,
                "No metrics rows available",
            ),
        ),
    ],
)
def test_build_tables(inputs, exception, wants, caplog):
    with exception:
        got = _build_tables(*inputs)
        for key in got.keys():
            if isinstance(wants[0][key], int):
                assert got[key] == wants[0][key]
            else:
                print(key)
                print(got[key])
                print(wants[0][key])
                assert got[key].equals(wants[0][key])
    assert wants[1] in caplog.text


def test_render_plot(caplog):
    """Unsupported renderer types raise ValueError."""
    with raises(ValueError):
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


def test_render_cluster_density_scatter(
    tmp_path,
):
    """Cluster-density renderer creates a non-empty PDF page."""
    output = tmp_path / "cluster_density.pdf"

    with PdfPages(output) as pdf_handle:
        _render_cluster_density_scatter(
            pdf_handle,
            {
                "source": "joint_qc_table",
                "title": "Cluster density",
                "skip_if_empty": True,
            },
            {
                "joint_qc_table": test_joint_qc_frame_table,
            },
            "dragen",
        )

    assert output.exists()
    assert output.stat().st_size > 1000
    assert output.read_bytes().startswith(b"%PDF")


def test_render_contamination_scatter(
    tmp_path,
):
    """Contamination renderer creates a non-empty PDF page."""
    output = tmp_path / "contamination.pdf"

    with PdfPages(output) as pdf_handle:
        _render_contamination_scatter(
            pdf_handle,
            {
                "source": "dna_data_table",
                "color_var": "RUN",
                "label_var": None,
                "title": "Contamination",
            },
            {
                "dna_data_table": pl.concat(
                    [
                        test_metrics_frame_table.slice(0, 1),
                        test_metrics_frame_table.slice(2, 1),
                    ]
                )
                .insert_column(9, pl.Series("highlighted_run", ["True", "False"]))
                .insert_column(
                    10, pl.Series("contamination_label", ["DNA_LATEST", ""])
                ),
                "dna_guideline_table": test_metrics_frame_guideline_table.slice(0, 2),
            },
            "localapp",
        )

    assert output.exists()
    assert output.stat().st_size > 1000
    assert output.read_bytes().startswith(b"%PDF")


@mark.parametrize(
    "inputs, exception, want",
    [
        (
            (
                test_metrics_frame_table.slice(0, 1),
                test_joint_qc_frame_table.slice(0, 1),
                "dragen",
            ),
            nullcontext(),
            "",
        ),
        (
            (
                test_metrics_frame_table.slice(0, 1).replace_column(
                    3,
                    pl.Series(
                        "WORKFLOW_TYPE",
                        ["localapp"],
                    ),
                ),
                test_joint_qc_frame_table.slice(0, 1).replace_column(
                    2,
                    pl.Series(
                        "WORKFLOW_TYPE",
                        ["localapp"],
                    ),
                ),
                "localapp",
            ),
            nullcontext(),
            "",
        ),
        (
            (
                test_metrics_frame_table,
                test_joint_qc_frame_table,
                "unknown",
            ),
            raises(ValueError),
            "Unsupported workflow",
        ),
    ],
)
def test_generate_qc_plots(
    inputs,
    exception,
    want,
    tmp_path,
    monkeypatch,
    caplog,
):
    """Unknown workflow values are rejected before rendering."""
    with exception:
        output_pdf = tmp_path / f"{input}_metric_plots.pdf"
        monkeypatch.setattr(
            plotting,
            "PLOT_SPECS",
            test_bar_plot_spec_complete,
        )
        Generate_qc_plots(
            metrics_table=inputs[0],
            joint_qc_table=inputs[1],
            workflow=inputs[2],
            output_pdf=output_pdf,
        )
        assert output_pdf.exists()
        assert output_pdf.stat().st_size > 1000

        content = output_pdf.read_bytes()

        assert content.startswith(b"%PDF")
        assert b"%%EOF" in content[-1024:]

    assert want in caplog.text
