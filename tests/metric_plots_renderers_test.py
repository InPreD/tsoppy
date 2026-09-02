"""Tests for the workflow-specific plotting renderers."""

from unittest.mock import MagicMock

import polars as pl
from matplotlib.backends.backend_pdf import PdfPages

import tsoppy.metric_plots.plotting as plotting
from tsoppy.metric_plots.plotting import (
    render_cluster_density_scatter,
    render_contamination_scatter,
)


def _cluster_data() -> pl.DataFrame:
    """Return representative sequencing run data."""
    return pl.DataFrame(
        {
            "RUN_ID": [
                "RUN_A",
                "RUN_B",
            ],
            "RUN_INDEX": [
                "001",
                "002",
            ],
            "CLUSTER_DENSITY": [
                "200",
                "250",
            ],
            "ESTIMATED_YIELD": [
                "100",
                "120",
            ],
        }
    )


def _contamination_data(
    max_score: str = "2000",
) -> pl.DataFrame:
    """Return representative DNA contamination data."""
    return pl.DataFrame(
        {
            "SAMPLE_ID": [
                "SAMPLE_A",
                "SAMPLE_B",
            ],
            "RUN": [
                "RUN_A",
                "RUN_B",
            ],
            "RUN_INDEX": [
                "001",
                "002",
            ],
            "DNA_CONTAMINATION_SCORE": [
                "100",
                max_score,
            ],
            "DNA_CONTAMINATION_P_VALUE": [
                "0.01",
                "0.10",
            ],
            "highlighted_run": [
                "True",
                "False",
            ],
            "contamination_label": [
                "SAMPLE_A",
                "",
            ],
        }
    )


def _dna_guideline_table(
    score: str = "1457",
    pvalue: str = "0.05",
) -> pl.DataFrame:
    """Return a representative DNA USL guideline."""
    return pl.DataFrame(
        {
            "SAMPLE_ID": [
                "USL_Guideline",
            ],
            "DNA_CONTAMINATION_SCORE": [
                score,
            ],
            "DNA_CONTAMINATION_P_VALUE": [
                pvalue,
            ],
        }
    )


def test_cluster_density_renderer_skips_empty_filtered_data(
    monkeypatch,
):
    """Cluster-density plots are skipped when no valid yield exists."""
    table = pl.DataFrame(
        {
            "RUN_ID": [
                "RUN_A",
                "RUN_B",
            ],
            "RUN_INDEX": [
                "001",
                "002",
            ],
            "CLUSTER_DENSITY": [
                "200",
                "250",
            ],
            "ESTIMATED_YIELD": [
                "NA",
                None,
            ],
        }
    )

    save_mock = MagicMock()

    monkeypatch.setattr(
        plotting,
        "save_plot",
        save_mock,
    )

    render_cluster_density_scatter(
        MagicMock(),
        {
            "source": "joint_qc_table",
            "title": "Cluster density",
            "skip_if_empty": True,
        },
        {
            "joint_qc_table": table,
        },
        "dragen",
    )

    save_mock.assert_not_called()


def test_cluster_density_renderer_creates_real_pdf(
    tmp_path,
):
    """Cluster-density renderer creates a non-empty PDF page."""
    output = tmp_path / "cluster_density.pdf"

    with PdfPages(output) as pdf_handle:
        render_cluster_density_scatter(
            pdf_handle,
            {
                "source": "joint_qc_table",
                "title": "Cluster density",
                "skip_if_empty": True,
            },
            {
                "joint_qc_table": _cluster_data(),
            },
            "dragen",
        )

    assert output.exists()
    assert output.stat().st_size > 1000
    assert output.read_bytes().startswith(b"%PDF")


def test_contamination_renderer_skips_empty_filtered_data(
    monkeypatch,
):
    """Contamination scatter is skipped when all values are missing."""
    table = pl.DataFrame(
        {
            "RUN": [
                "RUN_A",
            ],
            "RUN_INDEX": [
                "001",
            ],
            "DNA_CONTAMINATION_SCORE": [
                "NA",
            ],
            "DNA_CONTAMINATION_P_VALUE": [
                "NA",
            ],
        }
    )

    save_mock = MagicMock()

    monkeypatch.setattr(
        plotting,
        "save_plot",
        save_mock,
    )

    render_contamination_scatter(
        MagicMock(),
        {
            "source": "dna_data_table",
            "color_var": "RUN",
            "label_var": None,
            "title": "Contamination",
            "skip_if_empty": True,
        },
        {
            "dna_data_table": table,
        },
        "localapp",
    )

    save_mock.assert_not_called()


def test_contamination_renderer_uses_guideline_values(
    monkeypatch,
):
    """Configured contamination USLs are passed to the plot helper."""
    plot_object = MagicMock()

    plot_mock = MagicMock(
        return_value=plot_object,
    )

    save_mock = MagicMock()

    monkeypatch.setattr(
        plotting,
        "plot_contamination_scatter",
        plot_mock,
    )
    monkeypatch.setattr(
        plotting,
        "save_plot",
        save_mock,
    )

    render_contamination_scatter(
        MagicMock(),
        {
            "source": "dna_data_table",
            "color_var": "RUN",
            "label_var": None,
            "title": "Contamination",
        },
        {
            "dna_data_table": _contamination_data(),
            "dna_guideline_table": _dna_guideline_table(),
        },
        "localapp",
    )

    kwargs = plot_mock.call_args.kwargs

    assert kwargs["color_var"] == "RUN_LABEL"
    assert kwargs["usl_contamination_score"] == 1457.0
    assert kwargs["usl_contamination_pval"] == 0.05
    assert kwargs["max_contamination_score"] == 5000

    save_mock.assert_called_once()


def test_contamination_renderer_uses_default_values_for_na_guidelines(
    monkeypatch,
):
    """NA contamination guidelines use fallback plotting limits."""
    plot_mock = MagicMock(
        return_value=MagicMock(),
    )

    monkeypatch.setattr(
        plotting,
        "plot_contamination_scatter",
        plot_mock,
    )
    monkeypatch.setattr(
        plotting,
        "save_plot",
        MagicMock(),
    )

    render_contamination_scatter(
        MagicMock(),
        {
            "source": "dna_data_table",
            "color_var": "RUN",
            "label_var": None,
            "title": "Contamination",
        },
        {
            "dna_data_table": _contamination_data(),
            "dna_guideline_table": _dna_guideline_table(
                score="NA",
                pvalue="NA",
            ),
        },
        "localapp",
    )

    kwargs = plot_mock.call_args.kwargs

    assert kwargs["usl_contamination_score"] == 5000.0

    assert kwargs["usl_contamination_pval"] == 0.05


def test_contamination_renderer_expands_score_axis_above_5000(
    monkeypatch,
):
    """Observed scores above 5000 expand the contamination x range."""
    plot_mock = MagicMock(
        return_value=MagicMock(),
    )

    monkeypatch.setattr(
        plotting,
        "plot_contamination_scatter",
        plot_mock,
    )
    monkeypatch.setattr(
        plotting,
        "save_plot",
        MagicMock(),
    )

    render_contamination_scatter(
        MagicMock(),
        {
            "source": "dna_data_table",
            "color_var": "RUN",
            "label_var": None,
            "title": "Contamination",
        },
        {
            "dna_data_table": _contamination_data(
                max_score="7200",
            ),
            "dna_guideline_table": _dna_guideline_table(),
        },
        "localapp",
    )

    assert plot_mock.call_args.kwargs["max_contamination_score"] == 7200.0


def test_contamination_renderer_preserves_non_run_color_variable(
    monkeypatch,
):
    """Non-RUN color variables are passed through unchanged."""
    plot_mock = MagicMock(
        return_value=MagicMock(),
    )

    monkeypatch.setattr(
        plotting,
        "plot_contamination_scatter",
        plot_mock,
    )
    monkeypatch.setattr(
        plotting,
        "save_plot",
        MagicMock(),
    )

    render_contamination_scatter(
        MagicMock(),
        {
            "source": "dna_data_table",
            "color_var": "highlighted_run",
            "label_var": "contamination_label",
            "title": "Highlighted contamination",
            "color_values": {
                "False": "#888686",
                "True": "#C02F2F",
            },
        },
        {
            "dna_data_table": _contamination_data(),
            "dna_guideline_table": _dna_guideline_table(),
        },
        "localapp",
    )

    kwargs = plot_mock.call_args.kwargs

    assert kwargs["color_var"] == "highlighted_run"

    assert kwargs["label_var"] == "contamination_label"

    assert kwargs["color_values"] == {
        "False": "#888686",
        "True": "#C02F2F",
    }


def test_contamination_renderer_creates_real_pdf(
    tmp_path,
):
    """Contamination renderer creates a non-empty PDF page."""
    output = tmp_path / "contamination.pdf"

    with PdfPages(output) as pdf_handle:
        render_contamination_scatter(
            pdf_handle,
            {
                "source": "dna_data_table",
                "color_var": "RUN",
                "label_var": None,
                "title": "Contamination",
            },
            {
                "dna_data_table": _contamination_data(),
                "dna_guideline_table": _dna_guideline_table(),
            },
            "localapp",
        )

    assert output.exists()
    assert output.stat().st_size > 1000
    assert output.read_bytes().startswith(b"%PDF")
