"""
This module contains the CLI commands for tsoppy.
"""

import importlib.metadata
import logging
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


logger = logging.getLogger(__name__)


def validate_not_empty(file: Path) -> Path:
    """Callback function to ensure the file contains data."""
    # Typer runs built-in validation (like exists=True) BEFORE callbacks,
    # so we know the path already exists and is a file here.
    if file.stat().st_size == 0:
        logger.exception(f"The specified file {file} is empty.")
        raise typer.BadParameter()
    return file


@app.command()
def report_predispositions(
    output: Annotated[
        Path, typer.Option(help="File containing predisposition variants.")
    ],
    sample_id_tumor_DNA: Annotated[str, typer.Option(help="ID of the input tumor DNA sample.")],
    sample_id_normal_DNA: Annotated[str | None, typer.Option(
        help="ID of the input normal DNA sample.")] = None,
    sample_id_tumor_RNA: Annotated[str | None, typer.Option(
        help="ID of the input tumor RNA sample.")] = None,
    config_yaml: Annotated[Path, typer.Option(help="Config yaml file.")],
    inpred_nomenclature: Annotated[
        Path,
        typer.Option(
            exists=True,
            help="InPreD nomenclature definition file."
        )
    ],
    root_path_tumor_DNA: Annotated[
        Path,
        typer.Option(
            exists=True,
            help="Root path to the folder containing TSO500 analysis of the sequencing run which contains the input tumor DNA sample."
        ),
    ],
    root_path_normal_DNA: Annotated[] = None,
    root_path_tumor_RNA: Annotated[] = None,
    tumor_purity: Annotated[
        float,
        typer.Option(
            help="Proportion of cancerous cells in the analysed sample, the value has to be in the range between 0 and 1, where 1 represents 100%."
        ),
    ],
    cancer_susceptibility_genes: Annotated[
        Path,
        typer.Option(
            exists=True,
            callback=validate_not_empty,
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
        version(),
        sample_id_tumor_DNA,
        sample_id_normal_DNA,
        sample_id_tumor_RNA,
        config_yaml,
        inpred_nomenclature,
        root_path_tumor_DNA,
        root_path_normal_DNA,
        root_path_tumor_RNA,
        tumor_purity,
        cancer_susceptibility_genes,
        cancer_susceptibility_genes_column_list,
        cancer_susceptibility_genes_gene_name_column,
    )
