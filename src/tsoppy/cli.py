"""
This module contains the CLI commands for tsoppy.
"""

import importlib.metadata
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Annotated
from urllib import request
from tsoppy.general.classes import WorkflowOutput

import typer

app = typer.Typer()
app_version = importlib.metadata.version("tsoppy")


@dataclass
class CommonOptions:
    """Holds common input options."""

    nomenclature: Path
    config: Path


# TODO: Update url when PR is merged to main branch and uncomment this section
# def nomenclature_callback(value: Path | None) -> Path:
#    """
#    Callback function to provide fallback to nomenclature.yaml available on GitHub in case of missing or non-existent file.
#    """
#    if value is None or not value.exists():
#        temp_file = tempfile.NamedTemporaryFile(delete=False)
#        request.urlretrieve(
#            "https://raw.githubusercontent.com/InPreD/reference/refs/heads/5-add-inpred-nomenclature/InPreD/sample_id_nomenclature/v4/nomenclature.yaml",
#            temp_file.name,
#        )
#        return Path(temp_file.name)
#    else:
#        return value


@app.callback()
def main(
    ctx: typer.Context,
    config: Annotated[
        Path | None, typer.Option("--config", help="Path to tsoppy config file.")
    ] = None,
    nomenclature: Annotated[
        Path | None,
        typer.Option(
            "--nomenclature",
            help="Path to inpred nomenclature file.",
            # callback=nomenclature_callback,
        ),
    ] = None,
):
    """
    This is the main entry point for the tsoppy CLI application. It sets up the
    context object with commonOptions that can be accessed by subcommands.
    """
    ctx.obj = CommonOptions(nomenclature=nomenclature, config=config)


@app.command()
def version():
    """
    Prints the version of tsoppy.
    """
    print(f"tsoppy version {app_version}")


@app.command()
def placeholder(
    ctx: typer.Context,
):
    """
    This command illustrates how to include the commonOptions in a subcommand.
    """
    commonOptions = ctx.obj
    print(commonOptions.config)
    print(commonOptions.nomenclature)
