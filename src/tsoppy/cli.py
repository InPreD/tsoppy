"""
This module contains the CLI commands for tsoppy.
"""

import importlib.metadata
from typing import Annotated

import typer

from .report_predispositions import report_predispositions

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
):
    """
    This function reports variants called by small variant caller that are present 
    in the cancer susceptibility genes. Susceptibility of a gene depends on a patients age group, 
    thus age input parameter.  
    """
    
    report_predispositions.report_predispositions(cancer_susceptibility_genes,small_variant_calls,age,output)
    
