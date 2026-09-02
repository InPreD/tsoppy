"""End-to-end smoke tests for QC PDF generation."""

from unittest.mock import MagicMock

import polars as pl
import pytest
from typer.testing import CliRunner

import tsoppy.cli as cli_module
import tsoppy.metric_plots.plotting as plotting
from tsoppy.cli import app
from tsoppy.metric_plots.plotting import (
    Generate_qc_plots,
)

runner = CliRunner()


def _metrics_frame(
    workflow: str,
) -> pl.DataFrame:
    """Return the minimum metrics frame required by _build_tables."""
    return pl.DataFrame(
        {
            "SAMPLE_ID": [
                "DNA_SAMPLE_A",
            ],
            "RUN": [
                "RUN_A",
            ],
            "RUN_INDEX": [
                "001",
            ],
            "WORKFLOW_TYPE": [
                workflow,
            ],
            "RECORD_TYPE": [
                "DNA_SAMPLE",
            ],
            "DNA_CONTAMINATION_SCORE": [
                "100",
            ],
            "RNA_MEDIAN_CV_GENE_500X": [
                None,
            ],
        },
        schema_overrides={
            "RNA_MEDIAN_CV_GENE_500X": (pl.String),
        },
    )


def _joint_qc_frame(
    workflow: str,
) -> pl.DataFrame:
    """Return a minimum run-level QC frame."""
    return pl.DataFrame(
        {
            "RUN_ID": [
                "RUN_A",
            ],
            "RUN_INDEX": [
                "001",
            ],
            "WORKFLOW_TYPE": [
                workflow,
            ],
            "PCT_PF_READS": [
                "95",
            ],
        }
    )


def _minimal_plot_specs() -> dict:
    """Return one real bar plot enabled for both workflows."""
    return {
        "TEST_RUN_METRIC": {
            "localapp": {
                "plot": True,
                "index": 1,
            },
            "dragen": {
                "plot": True,
                "index": 1,
            },
            "plot_kind": "bar",
            "source": "joint_qc_table",
            "requires_samples": None,
            "title": "End-to-end QC test",
            "x_var": "RUN_ID",
            "y_var": "PCT_PF_READS",
            "fill_var": "RUN_ID",
            "x_lab": "Run ID",
            "y_lab": "Percentage",
            "value_spec": {
                "operation": "cast",
                "column": "PCT_PF_READS",
                "dtype": pl.Float64,
            },
            "na_filter_columns": [
                "PCT_PF_READS",
            ],
            "cart_ylim": (
                0,
                100,
            ),
            "skip_if_empty": True,
        }
    }


@pytest.mark.parametrize(
    "workflow",
    [
        "dragen",
        "localapp",
    ],
)
def test_generate_qc_plots_creates_real_pdf(
    monkeypatch,
    tmp_path,
    workflow,
):
    """Full plotting orchestration produces a real PDF."""
    monkeypatch.setattr(
        plotting,
        "PLOT_SPECS",
        _minimal_plot_specs(),
    )

    output = tmp_path / f"{workflow}_metric_plots.pdf"

    Generate_qc_plots(
        metrics_table=_metrics_frame(workflow),
        joint_qc_table=_joint_qc_frame(workflow),
        workflow=workflow,
        output_pdf=output,
    )

    assert output.exists()
    assert output.stat().st_size > 1000

    content = output.read_bytes()

    assert content.startswith(b"%PDF")
    assert b"%%EOF" in content[-1024:]


@pytest.mark.parametrize(
    "workflow",
    [
        "dragen",
        "localapp",
    ],
)
def test_cli_to_real_pdf_end_to_end(
    monkeypatch,
    tmp_path,
    workflow,
):
    """CLI selection flows through the real plotting stack to a PDF."""
    monkeypatch.chdir(tmp_path)

    config = tmp_path / "config.yaml"
    nomenclature = tmp_path / "nomenclature.yaml"

    config.write_text("{}\n")
    nomenclature.write_text("{}\n")

    metrics = _metrics_frame(workflow)
    joint_qc = _joint_qc_frame(workflow)

    metric_plotter = MagicMock()

    metric_plotter.generate_metrics_tables.return_value = (
        metrics,
        joint_qc,
    )

    metric_plotter.select_plot_data.return_value = (
        metrics,
        joint_qc,
    )

    monkeypatch.setattr(
        cli_module,
        "MetricPlots",
        MagicMock(
            return_value=metric_plotter,
        ),
    )

    monkeypatch.setattr(
        plotting,
        "PLOT_SPECS",
        _minimal_plot_specs(),
    )

    result = runner.invoke(
        app,
        [
            "metric-plots",
            "--input-glob",
            "results/*/*",
            "--config-yaml",
            str(config),
            "--inpred-nomenclature",
            str(nomenclature),
            "--run-ids",
            "RUN_A",
            "--plot-last-runs",
            "1",
            "--plot-workflow",
            workflow,
        ],
    )

    assert result.exit_code == 0, result.output

    metric_plotter.generate_metrics_tables.assert_called_once()

    metric_plotter.select_plot_data.assert_called_once()

    output = tmp_path / f"{workflow}_metric_plots.pdf"

    assert output.exists()
    assert output.stat().st_size > 1000

    content = output.read_bytes()

    assert content.startswith(b"%PDF")
    assert b"%%EOF" in content[-1024:]
