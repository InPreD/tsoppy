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
            "--input-directory",
            "-i",
            help="Input directory containing dragen/ and/or localapp/ MetricsOutput.tsv files.",
        ),
    ] = Path("in"),
    run_id_file: Annotated[
        Path,
        typer.Option(
            "--run-id-file",
            "-r",
            help="Text file containing one RUN ID per line.",
        ),
    ] = Path("in/RUN_IDs.txt"),
    output_directory: Annotated[
        Path,
        typer.Option(
            "--output-directory",
            "-o",
            help="Output directory for master table, joint QC file, and plots.",
        ),
    ] = Path("out"),
    create_plots: Annotated[
        bool,
        typer.Option(
            "--create-plots/--no-create-plots",
            help="Create PDF plots after generating metrics tables.",
        ),
    ] = True,
):
    """
    Create TSO500 metrics master table, joint sequencing QC file, and optional plots.
    """
    logger.info("Start metric plotting workflow.")

    metric_plotter = MetricPlots(
        input_directory=input_directory,
        run_id_file=run_id_file,
        output_directory=output_directory,
        create_plots=create_plots,
    )

    metric_plotter.run()

    logger.info("Finished metric plotting workflow.")
