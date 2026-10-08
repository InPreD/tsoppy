"""CLI module unit tests for tsoppy."""

from pathlib import Path
from unittest.mock import MagicMock

import polars as pl
import pytest
from click import unstyle
from typer.testing import CliRunner

import tsoppy.cli as cli_module
from tsoppy.cli import app

runner = CliRunner()


def test_placeholder():
    """Unit test for the placeholder command in the CLI module."""
    assert True


def _clean_output(result) -> str:
    """Remove terminal styling from CLI output."""
    return unstyle(result.output)


def _required_files(
    tmp_path,
) -> tuple[Path, Path]:
    """Create files required by Typer path validation."""
    config = tmp_path / "config.yaml"
    nomenclature = tmp_path / "nomenclature.yaml"

    config.write_text("{}\n")
    nomenclature.write_text("{}\n")

    return config, nomenclature


def _base_args(
    tmp_path,
) -> list[str]:
    """Return common valid CLI arguments."""
    config, nomenclature = _required_files(tmp_path)

    return [
        "metric-plots",
        "--input-glob",
        "results/*/*",
        "--config-yaml",
        str(config),
        "--inpred-nomenclature",
        str(nomenclature),
    ]


def _mock_metric_plotter(
    monkeypatch,
):
    """Install a mocked MetricPlots instance."""
    master = pl.DataFrame(
        {
            "MASTER": [1],
        }
    )

    joint_qc = pl.DataFrame(
        {
            "JOINT": [1],
        }
    )

    plot_frame = pl.DataFrame(
        {
            "PLOT": [1],
        }
    )

    plot_joint_qc = pl.DataFrame(
        {
            "PLOT_JOINT": [1],
        }
    )

    metric_plotter = MagicMock()

    metric_plotter.generate_metrics_tables.return_value = (
        master,
        joint_qc,
    )

    metric_plotter.select_plot_data.return_value = (
        plot_frame,
        plot_joint_qc,
    )

    constructor = MagicMock(
        return_value=metric_plotter,
    )

    monkeypatch.setattr(
        cli_module,
        "MetricPlots",
        constructor,
    )

    return {
        "constructor": constructor,
        "metric_plotter": metric_plotter,
        "master": master,
        "joint_qc": joint_qc,
        "plot_frame": plot_frame,
        "plot_joint_qc": plot_joint_qc,
    }


def test_cli_omits_master_run_selector_defers_to_input_glob(
    monkeypatch,
    tmp_path,
):
    """Omitting both master run selectors is valid and defers to --input-glob."""
    mocks = _mock_metric_plotter(monkeypatch)

    result = runner.invoke(
        app,
        _base_args(tmp_path),
    )

    assert result.exit_code == 0, result.output

    constructor_kwargs = mocks["constructor"].call_args.kwargs

    assert constructor_kwargs["run_ids"] is None


@pytest.mark.parametrize(
    "run_ids_first",
    [True, False],
)
def test_cli_rejects_both_master_run_selectors(
    run_ids_first,
    tmp_path,
):
    """--run-ids and --run-id-file cannot be combined, regardless of argument order."""
    run_file = tmp_path / "runs.txt"
    run_file.write_text("RUN_A\n")

    ids_args = ["--run-ids", "RUN_A"]
    file_args = ["--run-id-file", str(run_file)]
    extra_args = ids_args + file_args if run_ids_first else file_args + ids_args

    result = runner.invoke(
        app,
        _base_args(tmp_path) + extra_args,
    )

    assert result.exit_code != 0
    assert "--run-id-file" in _clean_output(result)


def test_cli_tables_only_does_not_generate_pdf(
    monkeypatch,
    tmp_path,
):
    """Table-only execution does not invoke the plotting layer."""
    mocks = _mock_metric_plotter(monkeypatch)

    generate_mock = MagicMock()

    monkeypatch.setattr(
        cli_module,
        "Generate_qc_plots",
        generate_mock,
    )

    result = runner.invoke(
        app,
        _base_args(tmp_path)
        + [
            "--run-ids",
            "RUN_A,RUN_B",
        ],
    )

    assert result.exit_code == 0, result.output

    mocks["metric_plotter"].generate_metrics_tables.assert_called_once()

    mocks["metric_plotter"].select_plot_data.assert_not_called()

    generate_mock.assert_not_called()


def test_cli_plotting_requires_workflow(
    tmp_path,
):
    """Plot selection requires --plot-workflow."""
    result = runner.invoke(
        app,
        _base_args(tmp_path)
        + [
            "--run-ids",
            "RUN_A",
            "--plot-last-runs",
            "2",
        ],
    )

    assert result.exit_code != 0
    assert "--plot-workflow" in _clean_output(result)


@pytest.mark.parametrize(
    "last_runs_first",
    [True, False],
)
def test_cli_rejects_last_runs_with_explicit_plot_runs(
    last_runs_first,
    tmp_path,
):
    """--plot-last-runs cannot be combined with explicit plot IDs, regardless of order."""
    last_runs_args = ["--plot-last-runs", "2"]
    plot_run_ids_args = ["--plot-run-ids", "RUN_A"]
    extra_args = (
        last_runs_args + plot_run_ids_args
        if last_runs_first
        else plot_run_ids_args + last_runs_args
    )

    result = runner.invoke(
        app,
        _base_args(tmp_path)
        + ["--run-ids", "RUN_A,RUN_B"]
        + extra_args
        + ["--plot-workflow", "dragen"],
    )

    assert result.exit_code != 0
    assert "--plot-last-runs" in _clean_output(result)


@pytest.mark.parametrize(
    "plot_run_ids_first",
    [True, False],
)
def test_cli_rejects_both_explicit_plot_run_selectors(
    plot_run_ids_first,
    tmp_path,
):
    """Plot run IDs cannot come from both CLI string and file, regardless of order."""
    plot_file = tmp_path / "plot_runs.txt"
    plot_file.write_text("RUN_A\n")

    plot_run_ids_args = ["--plot-run-ids", "RUN_A"]
    plot_run_id_file_args = ["--plot-run-id-file", str(plot_file)]
    extra_args = (
        plot_run_ids_args + plot_run_id_file_args
        if plot_run_ids_first
        else plot_run_id_file_args + plot_run_ids_args
    )

    result = runner.invoke(
        app,
        _base_args(tmp_path)
        + ["--run-ids", "RUN_A,RUN_B"]
        + extra_args
        + ["--plot-workflow", "dragen"],
    )

    assert result.exit_code != 0
    assert "--plot-run-id-file" in _clean_output(result)


def test_cli_plot_last_runs_must_be_positive(
    tmp_path,
):
    """Typer rejects zero as --plot-last-runs."""
    result = runner.invoke(
        app,
        _base_args(tmp_path)
        + [
            "--run-ids",
            "RUN_A",
            "--plot-last-runs",
            "0",
            "--plot-workflow",
            "dragen",
        ],
    )

    assert result.exit_code != 0


def test_cli_parses_and_deduplicates_plot_run_ids(
    monkeypatch,
    tmp_path,
):
    """Explicit plot-run strings are stripped and deduplicated."""
    mocks = _mock_metric_plotter(monkeypatch)

    generate_mock = MagicMock()

    monkeypatch.setattr(
        cli_module,
        "Generate_qc_plots",
        generate_mock,
    )

    result = runner.invoke(
        app,
        _base_args(tmp_path)
        + [
            "--run-ids",
            "RUN_A,RUN_B",
            "--plot-run-ids",
            " RUN_B, RUN_A,RUN_B ",
            "--plot-workflow",
            "dragen",
        ],
    )

    assert result.exit_code == 0, result.output

    kwargs = mocks["metric_plotter"].select_plot_data.call_args.kwargs

    assert kwargs["plot_run_ids"] == [
        "RUN_B",
        "RUN_A",
    ]

    assert kwargs["workflow_type"] == "dragen"


def test_cli_reads_plot_run_id_file(
    monkeypatch,
    tmp_path,
):
    """Plot-run files ignore blank lines and deduplicate IDs."""
    mocks = _mock_metric_plotter(monkeypatch)

    monkeypatch.setattr(
        cli_module,
        "Generate_qc_plots",
        MagicMock(),
    )

    plot_file = tmp_path / "plot_runs.txt"

    plot_file.write_text("\nRUN_B\nRUN_A\nRUN_B\n\n")

    result = runner.invoke(
        app,
        _base_args(tmp_path)
        + [
            "--run-ids",
            "RUN_A,RUN_B",
            "--plot-run-id-file",
            str(plot_file),
            "--plot-workflow",
            "localapp",
        ],
    )

    assert result.exit_code == 0, result.output

    kwargs = mocks["metric_plotter"].select_plot_data.call_args.kwargs

    assert kwargs["plot_run_ids"] == [
        "RUN_B",
        "RUN_A",
    ]

    assert kwargs["workflow_type"] == "localapp"


@pytest.mark.parametrize(
    ("workflow", "expected_filename"),
    [
        (
            "dragen",
            "dragen_metric_plots.pdf",
        ),
        (
            "localapp",
            "localapp_metric_plots.pdf",
        ),
    ],
)
def test_cli_passes_selected_frames_to_plot_generator(
    monkeypatch,
    tmp_path,
    workflow,
    expected_filename,
):
    """CLI passes selected frames and expected output path to plotting."""
    mocks = _mock_metric_plotter(monkeypatch)

    generate_mock = MagicMock()

    monkeypatch.setattr(
        cli_module,
        "Generate_qc_plots",
        generate_mock,
    )

    result = runner.invoke(
        app,
        _base_args(tmp_path)
        + [
            "--run-ids",
            "RUN_A,RUN_B",
            "--plot-last-runs",
            "2",
            "--plot-workflow",
            workflow,
        ],
    )

    assert result.exit_code == 0, result.output

    select_kwargs = mocks["metric_plotter"].select_plot_data.call_args.kwargs

    assert select_kwargs["plot_last_runs"] == 2

    assert select_kwargs["plot_run_ids"] is None

    assert select_kwargs["workflow_type"] == workflow

    generate_mock.assert_called_once()

    generate_kwargs = generate_mock.call_args.kwargs

    assert generate_kwargs["metrics_table"] is mocks["plot_frame"]

    assert generate_kwargs["joint_qc_table"] is mocks["plot_joint_qc"]

    assert generate_kwargs["workflow"] == workflow

    assert generate_kwargs["output_pdf"] == Path(expected_filename)
