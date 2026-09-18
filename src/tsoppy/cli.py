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


def _parse_run_ids(value: str) -> list[str]:
    """Split a comma-separated run ID string into stripped, non-empty IDs."""
    return [run_id.strip() for run_id in value.split(",") if run_id.strip()]


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
        # Annotated as str, not list[str], because typer/click treats a
        # list[str]-annotated option as repeatable and would wrap this
        # parser's own list return value in another list. The annotation
        # only controls CLI parsing arity; the parser is what actually
        # decides the value this function receives.
        str | None,
        typer.Option(
            parser=_parse_run_ids,
            metavar="TEXT",
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
            help=(
                "Text file containing run IDs for generation of the master metrics "
                "table, one per line. If neither --run-id-file nor --run-ids is "
                "provided, all runs matched by --input-glob are included. "
                "Mutually exclusive with --run-ids."
            ),
        ),
    ] = None,
    plot_run_ids: Annotated[
        # See the run_ids option above for why this stays annotated as str.
        str | None,
        typer.Option(
            parser=_parse_run_ids,
            metavar="TEXT",
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

    # The mutual-exclusivity checks below all run here in the function body
    # rather than as per-option callbacks: Click invokes option callbacks in
    # command-line argument order, not declaration order, so a callback
    # checking ctx.params for another option could silently miss a conflict
    # depending on the order options are typed on the command line.
    if run_ids is not None and run_id_file is not None:
        message = "--run-id-file cannot be used with --run-ids."
        logger.error(message)
        raise typer.BadParameter(message)

    if plot_run_ids is not None and plot_run_id_file is not None:
        message = "--plot-run-id-file cannot be used with --plot-run-ids."
        logger.error(message)
        raise typer.BadParameter(message)

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

    # run_ids is already parsed into a list by the option's parser=, if given.
    resolved_run_ids = run_ids

    if resolved_run_ids is None and run_id_file is not None:
        resolved_run_ids = [
            run_id
            for line in run_id_file.read_text().splitlines()
            if line.strip()
            for run_id in _parse_run_ids(line)
        ]

    metric_plotter = MetricPlots(
        config_yaml=config_yaml,
        inpred_nomenclature=inpred_nomenclature,
        input_glob=input_glob,
        run_ids=resolved_run_ids,
    )

    master, joint_qc = metric_plotter.generate_metrics_tables()

    logger.info("Metrics master table and joint QC files created.")

    if prepare_plot_frames:
        # plot_run_ids is already parsed into a list by the option's parser=, if given.
        plotting_run_ids = plot_run_ids

        if plotting_run_ids is None and plot_run_id_file is not None:
            plotting_run_ids = [
                line.strip()
                for line in plot_run_id_file.read_text().splitlines()
                if line.strip()
            ]

        if plotting_run_ids is not None:
            plotting_run_ids = list(dict.fromkeys(plotting_run_ids))

        if plot_workflow is None:
            logger.error(
                "Internal error: --plot-workflow was not resolved before "
                "plot selection despite passing the earlier validation."
            )
            raise RuntimeError

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
