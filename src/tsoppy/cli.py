"""
This module contains the CLI commands for tsoppy.
"""

import importlib.metadata
import logging
from typing import Annotated
from pathlib import Path

import typer

from tsoppy.metric_plots.main import MetricPlots


app = typer.Typer()

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s: %(message)s",
    datefmt="%Y/%m/%d %H:%M:%S",
)
logger = logging.getLogger(__name__)

app_version = importlib.metadata.version("tsoppy")


@app.command()
def version():
    """
    Prints the version of tsoppy.
    """
    print(f"tsoppy version {app_version}")


@app.command()
def placeholder(
    user_name: Annotated[str, typer.Option("--name", "-n")],
    user_id: Annotated[str, typer.Option("--id", "-i")],
    verbose: Annotated[bool, typer.Option("--verbose", "-v")] = False,
):
    """
    This is the helptext for the placeholder command that demonstrates how to
    use Typer for CLI applications.
    """
    if verbose:
        print(f"{user_name} has the following id: {user_id}")
    else:
        print(f"{user_name}: {user_id}")


@app.command()
def metric_plots(
    input_directory: Annotated[
        Path,
        typer.Option(
            help="Root directory containing workflow output directories.",
        ),
    ],
    config_yaml: Annotated[
        Path,
        typer.Option(
            help="Workflow configuration YAML.",
        ),
    ] = Path("config.yaml"),
    inpred_nomenclature: Annotated[
        Path,
        typer.Option(
            help="InPreD nomenclature YAML.",
        ),
    ] = Path("resources/nomenclature.yaml"),
    run_ids: Annotated[
        str | None,
        typer.Option(
            help=("Comma-separated run IDs to include in the master metrics table."),
        ),
    ] = None,
    run_id_file: Annotated[
        Path | None,
        typer.Option(
            help=(
                "Text file containing run IDs to include "
                "in the master metrics table, one per line."
            ),
        ),
    ] = None,
    workdir: Annotated[
        Path,
        typer.Option(
            help="Working directory for generated output files.",
        ),
    ] = Path("."),
    plot_last_runs: Annotated[
        int | None,
        typer.Option(
            help=(
                "Select the last N runs from the workflow "
                "specified with --plot-workflow."
            ),
        ),
    ] = None,
    plot_run_ids: Annotated[
        str | None,
        typer.Option(
            help=(
                "Comma-separated run IDs to select for plotting. "
                "Requires --plot-workflow."
            ),
        ),
    ] = None,
    plot_run_id_file: Annotated[
        Path | None,
        typer.Option(
            help=(
                "Text file containing run IDs to select for plotting. "
                "Requires --plot-workflow."
            ),
        ),
    ] = None,
    plot_workflow: Annotated[
        str | None,
        typer.Option(
            help="Workflow to plot: 'dragen' or 'localapp'.",
        ),
    ] = None,
):
    """Create metrics tables and optionally prepare data for plotting."""

    logger.info("Start creating metrics tables.")

    # Master-table run IDs
    master_run_ids = []

    if run_ids:
        master_run_ids.extend(
            run_id.strip() for run_id in run_ids.split(",") if run_id.strip()
        )

    if run_id_file:
        master_run_ids.extend(
            line.strip()
            for line in run_id_file.read_text().splitlines()
            if line.strip() and not line.strip().startswith("#")
        )

    master_run_ids = list(dict.fromkeys(master_run_ids))

    if not master_run_ids:
        raise typer.BadParameter(
            "Provide run IDs for the master table using --run-ids or --run-id-file."
        )

    # Plot-selection mode
    explicit_plot_runs = plot_run_ids is not None or plot_run_id_file is not None

    plotting_requested = plot_last_runs is not None or explicit_plot_runs

    if plot_last_runs is not None and explicit_plot_runs:
        raise typer.BadParameter(
            "--plot-last-runs cannot be combined with "
            "--plot-run-ids or --plot-run-id-file."
        )

    if plotting_requested and plot_workflow is None:
        raise typer.BadParameter(
            "--plot-workflow is required when plotting is requested."
        )

    if plot_workflow is not None and plot_workflow not in {"dragen", "localapp"}:
        raise typer.BadParameter(
            "--plot-workflow must be either 'dragen' or 'localapp'."
        )

    if plot_last_runs is not None and plot_last_runs < 1:
        raise typer.BadParameter("--plot-last-runs must be greater than zero.")

    metric_plotter = MetricPlots(
        config_yaml=config_yaml,
        inpred_nomenclature=inpred_nomenclature,
        input_directory=input_directory,
        workdir=workdir,
        run_ids=run_ids,
        run_id_file=run_id_file,
    )

    master, joint_qc = metric_plotter.run()

    if plotting_requested:
        plotting_run_ids = []

        if plot_run_ids:
            plotting_run_ids.extend(
                run_id.strip() for run_id in plot_run_ids.split(",") if run_id.strip()
            )

        if plot_run_id_file:
            plotting_run_ids.extend(
                line.strip()
                for line in plot_run_id_file.read_text().splitlines()
                if line.strip() and not line.strip().startswith("#")
            )

        plotting_run_ids = list(dict.fromkeys(plotting_run_ids))

        plot_frame, plot_joint_qc = metric_plotter.prepare_plot_frames(
            master=master,
            joint_qc=joint_qc,
            workflow_type=plot_workflow,
            plot_last_runs=plot_last_runs,
            plot_run_ids=(plotting_run_ids if plotting_run_ids else None),
        )

        #        plot_frame.write_csv(
        #            "plot_frame.tsv",
        #            separator="\t",
        #        )

        logger.info("Plot frame:\n%s", plot_frame)

        # Future:
        # create_qc_plots(
        #     metrics=plot_frame,
        #    joint_qc=plot_joint_qc,
        #     workflow_type=plot_workflow,
        #     output="metric_plots.pdf",
        # )

        logger.info(
            "Prepared %d rows for %s plotting.",
            plot_frame.height,
            plot_workflow,
        )

    logger.info("Finished creating metrics tables.")
