"""
This module contains the CLI commands for tsoppy.
"""

import importlib.metadata
import logging
import re
from pathlib import Path
from typing import Annotated

import typer

from tsoppy.update_small_variant_vcf_list.main import VcfList

# Set up logging for the CLI. The logging level is set to INFO, and the log messages will include the timestamp, log level, and message.
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(levelname)s: %(message)s',
                    datefmt='%Y/%m/%d %H:%M:%S')
logger = logging.getLogger(__name__)

app = typer.Typer()
app_version = importlib.metadata.version("tsoppy")


@app.command()
def version():
    """
    Prints the version of tsoppy.
    """
    print(f"tsoppy version {app_version}")


def glob_pattern_callback(value: str) -> str:
    """
    Callback function checking that the glob pattern ends with '.vcf'.
    """
    if not value.endswith(".vcf"):
        raise typer.BadParameter("Glob pattern must end with '.vcf'.")
    return value


def inpred_id_regex_callback(value: str) -> str:
    """
    Callback function ensuring inpred_id_regex contains the required named capture groups.
    """
    if "<patient_id>" not in value:
        raise typer.BadParameter(
            "inpred_id_regex must contain a named group 'patient_id'.")
    if "<sample_type>" not in value:
        raise typer.BadParameter(
            "inpred_id_regex must contain a named group 'sample_type'.")
    return value


def tumor_sample_types_callback(value: str) -> str:
    """
    Callback function to ensure tumor_sample_types is a comma-separated list of single letters.
    """
    if not re.fullmatch(r"^([A-Za-z],)+[A-Za-z]$", value):
        raise typer.BadParameter(
            "tumor_sample_types must be comma-separated list of single letters.")
    return value


@app.command()
def update_small_variant_vcf_list(
    results_dir: Annotated[Path | None, typer.Option(help="Directory where the results of the latest TSO500 run are stored.")],
    glob_pattern: Annotated[str, typer.Option(
        help="Glob pattern to search for small variant VCF files in the results directory.", callback=glob_pattern_callback)] = "**/Results/**/*_MergedSmallVariants.genome.vcf",
    inpred_id_regex: Annotated[str, typer.Option(
        help="Regular expression to extract the inpred_id from the VCF file name.", callback=inpred_id_regex_callback)] = "(?P<patient_id>\D{3}\d{4})-\D\d{2}-(?P<sample_type>\D)\d{2}-\D\d{2}.*.vcf$",
    output: Annotated[str, typer.Option(
        help="Name of new small variant VCF list.")] = f"small_variant_vcf_list_<YYYYMMDD>.tsv",
    tumor_sample_types: Annotated[str, typer.Option(
        help="Comma-separated list of sample types that are considered tumor samples.")] = "C,D,d,L,M,P,p,R,r,T,X",
    vcf_list: Annotated[Path | None, typer.Option(
        help="Path to list of small variant VCF files.")] = None,
):
    """
    Updates the small variant VCF list based on VCF(s) in results directory.
    """
    logger.info("Start updating small variant VCF list.")
    small_variant_vcf_list = VcfList(
        results_dir, glob_pattern, vcf_list, inpred_id_regex, tumor_sample_types, output)
    small_variant_vcf_list.update()
    logger.info("Finished updating small variant VCF list.")
