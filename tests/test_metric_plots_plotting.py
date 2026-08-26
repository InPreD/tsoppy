"""Unit tests for workflow-specific metric plotting helpers."""

import polars as pl
import pytest

from tsoppy.metric_plots.plotting import (
    build_tables,
    prepare_bar_plot_data,
    resolve_plot_title,
    validate_plot_specs,
)
from tsoppy.metric_plots.specs.plot_specs_workflows import PLOT_SPECS


def _metrics_frame() -> pl.DataFrame:
    """Build synthetic metrics rows for plotting-table tests."""
    return pl.DataFrame(
        {
            "SAMPLE_ID": [
                "DNA_LATEST",
                "RNA_LATEST",
                "DNA_OLDER",
                "LSL_Guideline",
                "USL_Guideline",
            ],
            "RUN": [
                "RUN_002",
                "RUN_002",
                "RUN_010",
                "GUIDELINE",
                "GUIDELINE",
            ],
            "RUN_INDEX": ["002", "002", "010", "999", "999"],
            "WORKFLOW_TYPE": ["dragen"] * 5,
            "RECORD_TYPE": [
                "DNA_SAMPLE",
                "RNA_SAMPLE",
                "DNA_SAMPLE",
                "LOWER_THRESHOLD",
                "UPPER_THRESHOLD",
            ],
            "DNA_CONTAMINATION_SCORE": [
                "100",
                "200",
                "110",
                "0",
                "1457",
            ],
            "RNA_MEDIAN_CV_GENE_500X": [
                "0.3",
                "0.2",
                "0.4",
                "0.0",
                "1.0",
            ],
        }
    )


def _joint_qc_frame() -> pl.DataFrame:
    """Build synthetic joint-QC rows for plotting-table tests."""
    return pl.DataFrame(
        {
            "RUN_ID": [
                "RUN_002",
                "RUN_010",
                "LSL_Guideline",
            ],
            "RUN_INDEX": [
                "002",
                "010",
                "999",
            ],
            "WORKFLOW_TYPE": ["dragen"] * 3,
        }
    )


def _minimal_bar_spec(index: int) -> dict:
    """Return a structurally valid minimal bar specification."""
    return {
        "localapp": {"plot": True, "index": index},
        "dragen": {"plot": False, "index": 0},
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


def test_run_index_sample_id_generation_and_padding():
    """Sample labels include RUN_INDEX and are padded to equal width."""
    table = pl.DataFrame(
        {
            "RUN_INDEX": ["001", "002"],
            "SAMPLE_ID": ["S1", "LONG_SAMPLE"],
            "RUN": ["RUN_A", "RUN_B"],
            "VALUE": [1.0, 2.0],
        }
    )

    plot_data = prepare_bar_plot_data(
        table,
        _minimal_bar_spec(index=1),
    )

    labels = plot_data["PLOT_SAMPLE_ID"].tolist()

    assert labels[0].rstrip() == "001 | S1"
    assert labels[1].rstrip() == "002 | LONG_SAMPLE"
    assert len(labels[0]) == len(labels[1])


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

    plot_data = prepare_bar_plot_data(table, spec)

    assert plot_data["PLOT_RUN"].tolist() == [
        "001 | RUN_A",
        "002 | RUN_B",
    ]


def test_latest_run_uses_minimum_numeric_run_index():
    """The latest highlighted run is the minimum numeric RUN_INDEX."""
    tables = build_tables(
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


def test_record_type_controls_dna_rna_selection():
    """DNA/RNA plotting tables are selected from RECORD_TYPE."""
    tables = build_tables(
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


def test_current_plot_specs_validate():
    """The production plot specification set is structurally valid."""
    validate_plot_specs(PLOT_SPECS)


def test_duplicate_plot_indices_are_rejected():
    """Duplicate positive plot indices within one workflow are rejected."""
    duplicate_specs = {
        "FIRST": _minimal_bar_spec(index=1),
        "SECOND": _minimal_bar_spec(index=1),
    }

    with pytest.raises(
        KeyError,
        match="duplicate plot index 1",
    ):
        validate_plot_specs(duplicate_specs)


def test_dragen_contamination_scatter_is_disabled():
    """LocalApp contamination scatter plots remain disabled for DRAGEN."""
    for spec_name in (
        "DNA_CONTAMINATION_P_VALUE_BY_RUN",
        "DNA_CONTAMINATION_P_VALUE_HIGHLIGHTED_RUN",
    ):
        assert PLOT_SPECS[spec_name]["dragen"] == {
            "plot": False,
            "index": 0,
        }


def test_dragen_contamination_score_bar_is_enabled():
    """DRAGEN uses the contamination-score bar plot at index 10."""
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
    """The five run-level DRAGEN plots remain enabled at indices 1-5."""
    expected = {
        "CLUSTERS_PASSING_FILTER": 1,
        "ESTIMATED_YIELD": 2,
        "PCT_PF_READS": 3,
        "PCT_Q30_R1": 4,
        "PCT_Q30_R2": 5,
    }

    for spec_name, expected_index in expected.items():
        assert PLOT_SPECS[spec_name]["dragen"] == {
            "plot": True,
            "index": expected_index,
        }


def test_workflow_specific_title_resolution():
    """Workflow-specific titles resolve correctly."""
    q30_spec = PLOT_SPECS["PCT_Q30_R1"]

    assert resolve_plot_title(
        q30_spec,
        "localapp",
    ).startswith("[LocalApp")

    assert resolve_plot_title(
        q30_spec,
        "dragen",
    ).startswith("[Dragen")

    shared_spec = PLOT_SPECS["CLUSTERS_PASSING_FILTER"]

    assert (
        resolve_plot_title(
            shared_spec,
            "localapp",
        )
        == shared_spec["title"]
    )

    assert (
        resolve_plot_title(
            shared_spec,
            "dragen",
        )
        == shared_spec["title"]
    )
