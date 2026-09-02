"""This module contains the CLI commands for tsoppy."""

import importlib.metadata
import logging
from enum import Enum
from pathlib import Path
from typing import Annotated

import typer

from tsoppy.metric_plots.main import MetricPlots
from tsoppy.metric_plots.plotting import Generate_qc_plots


class WorkflowType(str, Enum):
    """Supported workflow types for metric plotting."""

    DRAGEN = "dragen"
    LOCALAPP = "localapp"


app = typer.Typer()

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s: %(message)s",
    datefmt="%Y/%m/%d %H:%M:%S",
)
logger = logging.getLogger(__name__)

app_version = importlib.metadata.version("tsoppy")


def validate_run_id_file(
    ctx: typer.Context,
    value: Path | None,
) -> Path | None:
    """Ensure --run-id-file is not used with --run-ids."""
    if value is not None and ctx.params.get("run_ids") is not None:
        message = "--run-id-file cannot be used with --run-ids."
        logger.error(message)
        raise typer.BadParameter(message)
    return value


def validate_plot_run_id_file(
    ctx: typer.Context,
    value: Path | None,
) -> Path | None:
    """Ensure --plot-run-id-file is not used with --plot-run-ids."""
    if value is not None and ctx.params.get("plot_run_ids") is not None:
        message = "--plot-run-id-file cannot be used with --plot-run-ids."
        logger.error(message)
        raise typer.BadParameter(message)
    return value


@app.command()
def version():
    """Print the version of tsoppy."""
    print(f"tsoppy version {app_version}")


@app.command()
def placeholder(
    user_name: Annotated[
        str,
        typer.Option("--name", "-n"),
    ],
    user_id: Annotated[
        str,
        typer.Option("--id", "-i"),
    ],
    verbose: Annotated[
        bool,
        typer.Option("--verbose", "-v"),
    ] = False,
):
    """Demonstrate how to use Typer for CLI applications."""
    if verbose:
        print(f"{user_name} has the following id: {user_id}")
    else:
        print(f"{user_name}: {user_id}")


@app.command()
def metric_plots(
    input_glob: Annotated[
        str,
        typer.Option(
            help=(
                "Glob pattern matching workflow output directories "
                "whose final directory name is the sequencing run ID. "
                "Example:"
                " 'results/*/*' "
            ),
        ),
    ],
    inpred_nomenclature: Annotated[
        Path,
        typer.Option(
            exists=True,
            file_okay=True,
            dir_okay=False,
            readable=True,
            resolve_path=True,
            help="InPreD nomenclature YAML.",
        ),
    ],
    config_yaml: Annotated[
        Path,
        typer.Option(
            exists=True,
            file_okay=True,
            dir_okay=False,
            readable=True,
            resolve_path=True,
            help="Workflow configuration YAML.",
        ),
    ] = Path("config.yaml"),
    run_ids: Annotated[
        str | None,
        typer.Option(
            help=(
                "Comma-separated run IDs to include in the generated master metrics "
                "table. If neither --run-ids nor --run-id-file is provided, all runs "
                "matched by --input-glob are included. Mutually exclusive with "
                "--run-id-file."
            ),
        ),
    ] = None,
    run_id_file: Annotated[
        Path | None,
        typer.Option(
            exists=True,
            file_okay=True,
            dir_okay=False,
            readable=True,
            resolve_path=True,
            callback=validate_run_id_file,
            help=(
                "Text file containing run IDs for generation of the master metrics "
                "table, one per line. If neither --run-id-file nor --run-ids is "
                "provided, all runs matched by --input-glob are included. "
                "Mutually exclusive with --run-ids."
            ),
        ),
    ] = None,
    plot_run_ids: Annotated[
        str | None,
        typer.Option(
            help=(
                "Comma-separated list of run IDs to include in plot. "
                "Mutually exclusive with --plot-run-id-file."
            ),
        ),
    ] = None,
    plot_run_id_file: Annotated[
        Path | None,
        typer.Option(
            exists=True,
            file_okay=True,
            dir_okay=False,
            readable=True,
            resolve_path=True,
            callback=validate_plot_run_id_file,
            help=(
                "Text file containing list of run IDs to select for plotting. "
                "Mutually exclusive with --plot-run-ids."
            ),
        ),
    ] = None,
    plot_last_runs: Annotated[
        int | None,
        typer.Option(
            min=1,
            help=(
                "Plot the most recent N runs for the selected workflow. "
                "Mutually exclusive with --plot-run-ids and "
                "--plot-run-id-file."
            ),
        ),
    ] = None,
    plot_workflow: Annotated[
        WorkflowType | None,
        typer.Option(
            help=(
                "Workflow whose runs will be plotted. If no plot run selector is "
                "provided, the last 10 runs for this workflow are plotted."
            ),
        ),
    ] = None,
):
    """Create metrics tables and optionally generate QC plots."""
    logger.info("Creating metrics master table and joint QC.")

    plot_run_selection_given = plot_run_ids is not None or plot_run_id_file is not None

    prepare_plot_frames = (
        plot_workflow is not None
        or plot_last_runs is not None
        or plot_run_selection_given
    )

    # --plot-last-runs is a separate selection mode.
    if plot_last_runs is not None and plot_run_selection_given:
        message = (
            "--plot-last-runs cannot be combined with "
            "--plot-run-ids or --plot-run-id-file."
        )
        logger.error(message)
        raise typer.BadParameter(message)

    if prepare_plot_frames and plot_workflow is None:
        message = "--plot-workflow is required when plotting is requested."
        logger.error(message)
        raise typer.BadParameter(message)

    metric_plotter = MetricPlots(
        config_yaml=config_yaml,
        inpred_nomenclature=inpred_nomenclature,
        input_glob=input_glob,
        run_ids=run_ids,
        run_id_file=run_id_file,
    )

    master, joint_qc = metric_plotter.generate_metrics_tables()

    logger.info("Metrics master table and joint QC files created.")

    if prepare_plot_frames:
        plotting_run_ids: list[str] | None = None

        if plot_run_ids is not None:
            plotting_run_ids = [
                run_id.strip() for run_id in plot_run_ids.split(",") if run_id.strip()
            ]

        elif plot_run_id_file is not None:
            plotting_run_ids = [
                line.strip()
                for line in plot_run_id_file.read_text().splitlines()
                if line.strip() and not line.strip().startswith("#")
            ]

        if plotting_run_ids is not None:
            plotting_run_ids = list(dict.fromkeys(plotting_run_ids))

        if plot_workflow is None:
            message = (
                "Internal error: --plot-workflow was not resolved before "
                "plot selection despite passing the earlier validation."
            )
            logger.error(message)
            raise RuntimeError(message)

        plot_frame, plot_joint_qc = metric_plotter.select_plot_data(
            master=master,
            joint_qc=joint_qc,
            workflow_type=plot_workflow.value,
            plot_last_runs=plot_last_runs,
            plot_run_ids=plotting_run_ids,
        )

        logger.info(
            f"Prepared {plot_frame.height} metric rows and "
            f"{plot_joint_qc.height} joint QC rows for {plot_workflow.value} plotting."
        )

        output_pdf = Path(f"{plot_workflow.value}_metric_plots.pdf")

        Generate_qc_plots(
            metrics_table=plot_frame,
            joint_qc_table=plot_joint_qc,
            workflow=plot_workflow.value,
            output_pdf=output_pdf,
        )
