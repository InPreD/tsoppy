"""
This module contains the CLI commands for tsoppy.
"""

import importlib.metadata
from pathlib import Path
from typing import Annotated

import typer

from tsoppy.report_predispositions.main import generate_report

app = typer.Typer()
app_version = importlib.metadata.version("tsoppy")


@app.command()
def version():
    """
    Prints the version of tsoppy.
    """
    print(f"tsoppy version {app_version}")


@app.command()
def report_predispositions(
    output: Annotated[
        Path, typer.Option(help="File containing predisposition variants.")
    ],
    sample_id: Annotated[str, typer.Option(help="ID of the input sample.")],
    config_yaml: Annotated[Path, typer.Option(help="Config yaml file.")],
    inpred_nomenclature: Annotated[
        Path, typer.Option(help="InPreD nomenclature definition file.")
    ],
    root_path: Annotated[
        Path,
        typer.Option(
            help="Root path to the folder containing TSO500 analysis the sequencing run which contains the input sample."
        ),
    ],
    tumor_purity: Annotated[
        float,
        typer.Option(
            help="Proportion of cancerous cells in the analysed sample, the value has to be in the range between 0 and 1, where 1 represents 100%."
        ),
    ],
    cancer_susceptibility_genes: Annotated[
        Path,
        typer.Option(
            help="File containing info about genes in which variants increase a person's risk of developing certain cancers."
        ),
    ],
    cancer_susceptibility_genes_column_list: Annotated[
        list[str],
        typer.Option(
            help="List of column names for the columns that should be reported in the output predisposition file."
        ),
    ] = ["Gene", "Actionability", "Age"],
    cancer_susceptibility_genes_gene_name_column: Annotated[
        str,
        typer.Option(
            help="Name of the column containing HUGO gene names (e.g. BRCA1, TP53)."
        ),
    ] = "Gene",
):
    """
    This function reports variants called by small variant caller that are present
    in the cancer susceptibility genes.
    """

    generate_report(
        output,
        app_version,
        sample_id,
        config_yaml,
        inpred_nomenclature,
        root_path,
        tumor_purity,
        cancer_susceptibility_genes,
        cancer_susceptibility_genes_column_list,
        cancer_susceptibility_genes_gene_name_column,
    )
