"""
This module contains the CLI commands for tsoppy.
"""

import importlib.metadata
from typing import Annotated

import typer

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
def report_predispositions(
        cancer_susceptibility_genes: Annotated[str, typer.Option("--cancer-susceptibility-genes", "-g")],
        small_variant_calls: Annotated[str, typer.Option("--small-variant-calls", "-c")],
        age: Annotated[int, typer.Option("--age", "-a")],
        output: Annotated[str, typer.Option("--output", "-o") ],
        verbose: Annotated[bool, typer.Option("--verbose", "-v")] = False,
):
    """
    This function reports variants called by small variant caller that are present 
    in the cancer susceptibility genes.    
    """
    
    predispositions = dict()

    if verbose:
        print(f"Load data from the {cancer_susceptibility_genes} table")
    if verbose:
        print(f"Open the {small_variant_calls} file and iterate through the variants. 
                Store all the variants present in the {cancer_susceptibility_genes} table 
                together with all the info that should be reported into the {predispositions}.")
    if verbose:
        print(f"Print the {predispositions} content into the {output} file.")
