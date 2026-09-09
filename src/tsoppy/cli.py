"""
This module contains the CLI commands for tsoppy.
"""

import importlib.metadata
from pathlib import Path
from typing import Annotated

import typer

from tsoppy.add_variant_summary.main import Summarise_variant_output

app = typer.Typer()
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
def add_variant_summary(
    config_yaml: Annotated[
        Path, typer.Option("--config-yaml", "-c", help="Path to the config file.")
    ],
    inpred_nomenclature: Annotated[
        Path,
        typer.Option(
            "--inpred-nomenclature", "-i", help="Path to the inpred nomenclature file."
        ),
    ],
    results_directory: Annotated[
        Path,
        typer.Option(
            "--results-directory",
            "-r",
            help="path to the directory containing *_CombinedVariantOutput.tsv files.",
        ),
    ],
    output_file: Annotated[
        Path, typer.Option("--output-file", "-o", help="Path to the output file.")
    ],
):
    """Parse *_CombinedVariantOutput.tsv files into a variant summary table."""
    Summarise_variant_output(
        config_yaml, inpred_nomenclature, results_directory, output_file
    )
