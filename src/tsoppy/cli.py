"""
This module contains the CLI commands for tsoppy.
"""

import importlib.metadata
from pathlib import Path
from typing import Annotated

import typer

from tsoppy.report_predispositions.main import report_predispositions

app = typer.Typer()
app_version = importlib.metadata.version("tsoppy")


@app.command()
def version():
    """
    Prints the version of tsoppy.
    """
    print(f"tsoppy version {app_version}")


@app.command()
def report_predisposition_variants(
    sample_id: Annotated[str, typer.Option(help="ID of the input sample.")],
    version_string: Annotated[
        str,
        typer.Option(
            help="String of the tools (and their versions) that were used in the data analysis."
        ),
    ],
    doi_reference: Annotated[
        str,
        typer.Option(
            help="Url of the article recommanding the cancer susceptibility gene list that is defined in the '--cancer-susceptibility-genes' option."
        ),
    ],
    length_of_targeted_coding_regions: Annotated[
        float,
        typer.Option(
            help="Cummulative length of all coding regions targeted in the TSO500 experiment (in millions of bases)."
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
            help="List of column names for the columns that should be reported in the output predisposition file. The first column has to contain gene names (e.g. BRCA1, TP53)."
        ),
    ],
    germline_small_variant_calls: Annotated[
        Path, typer.Option(help="File containing germline small variant calls.")
    ],
    output: Annotated[
        Path, typer.Option(help="File containing predisposition variants.")
    ],
):
    """
    This function reports variants called by small variant caller that are present
    in the cancer susceptibility genes.
    """

    report_predispositions(
        sample_id,
        version_string,
        doi_reference,
        length_of_targeted_coding_regions,
        tumor_purity,
        cancer_susceptibility_genes,
        cancer_susceptibility_genes_column_list,
        germline_small_variant_calls,
        output,
    )
