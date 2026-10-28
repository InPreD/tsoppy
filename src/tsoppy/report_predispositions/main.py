"""
This module filters an input list of small variants and outputs a table containing only variants located in cancer-predisposing genes.
"""

import logging
import re
from collections import defaultdict
from pathlib import Path
from polars.exceptions import (ColumnNotFoundError, NoDataError)
from string import Template

import polars as pl
import sys

from tsoppy.general.classes import (
    SmallVariantGenomeVcf,
    WorkflowOutput
)

logger = logging.getLogger(__name__)


def validate_uniqueness(df: pl.DataFrame, column: str):
    """
    Test that all the values in column 'column' of dataframe 'df' are unique.
    """

    logger.info(
        f"Check that column \"{column}\" is a primary key of the data frame.")

    if not df[column].is_unique().all():
        # find the duplicates to make the error message helpful
        duplicates = df[column].filter(
            df[column].is_duplicated()).unique().to_list()
        # this captures the full traceback automatically
        logger.exception("Data Integrity Validation Failed.")
        raise ValueError(f"Duplicate IDs found in '{column}': {duplicates}")
    else:
        return


def validate_tumor_purity_range(tumor_purity: float):
    """
    Test that tumor purity value is between 0 and 1.
    """

    logger.info(
        f"Test that tumor purity value {tumor_purity} is between 0 and 1.")

    if not 0 <= tumor_purity <= 1:
        # this captures the full traceback automatically
        logger.exception("Data Integrity Validation Failed.")
        raise ValueError(f"Tumor purity value {tumor_purity} is out of range.")
    else:
        return


def check_gene_name_column_presence_in_column_list(gene_name_column: str, column_list: list[str]):
    """
    Check that gene_name_column is present in the provided column_list.
    """

    logger.info(
        f"Check that column \"{gene_name_column}\" is present in column list {column_list}.")

    if not gene_name_column in column_list:
        logger.exception("Data Integrity Validation Failed.")
        raise ValueError(
            f"Column {gene_name_column} is not present in provided column list {column_list}.")
    else:
        return


def check_file_existence_and_content_presence(file_path: Path):
    """
    Check that input file exists and is not empty.
    """

    logger.info("Check that input file {file_path} exists and is not empty.")

    if not file_path.is_file():
        logger.error("File Not Found.")
        raise FileNotFoundError(
            f"Error: file that should contain cancer susceptibility genes info ({file_path}) does not exist.")

    if file_path.stat().st_size == 0:
        logger.error("Empty Input File.")
        raise NoDataError(f"Error: Input file {file_path} is empty.")

    return


def parse_metadata_to_get_source(file_path: Path
                                 ) -> str:
    """
    Parse metadata of the file to get source of the information stored in the file.
    """

    logger.info(
        f"Parsing metadata of {file_path} to get source of the information stored in the file.")

    # try to get metadata about source of the info
    source = "NA"

    with open(file_path, "r") as file:
        for line in file:
            if re.match(r"# source:", line):
                source = line.replace("# source:", "", 1)
                source = source.strip()
                return source

    logger.exception("Data Integrity Validation Failed.")
    raise ValueError(
        f"Error: source information not found in metadata of {file_path}.")


def check_columns_presence_in_file(column_list: list[str], file_path: Path):
    """
    Check that all requested columns are present in the input file.
    """

    logger.info(
        f"Check that all requested columns from \"{column_list}\" are present in file {file_path}.")

    lf = pl.scan_csv(file_path, separator="\t", comment_prefix="#")

    available_columns = lf.columns
    missing_columns = [
        col for col in column_list if col not in available_columns]

    if missing_columns:
        logger.error(
            "Failed to read CSV because a required column was missing.")
        raise ColumnNotFoundError(f"Error: missing columns {missing_columns}.")
    else:
        return


def load_data_from_cancer_susceptibility_genes_file(
    file_path: Path, column_list: list[str], gene_name_column: str
) -> tuple[dict[str, dict[str, str]], str]:
    """
    Load data from the input cancer susceptibility genes file.
    """

    logger.info(
        f"Load data from the cancer susceptibility genes file {file_path}.")

    # check that gene_name_column is present in the column_list
    check_gene_name_column_presence_in_column_list(
        gene_name_column, column_list)

    # check the file existence and content presence
    check_file_existence_and_content_presence(file_path)

    # get intormation about the source of data presented in the input file
    source = parse_metadata_to_get_source(file_path)

    # check that all columns from the column list are present in the file
    check_columns_presence_in_file(column_list, file_path)

    logger.info(f"Read tabular content of the file {file_path}.")
    df = pl.read_csv(file_path, columns=column_list,
                     separator="\t", comment_prefix="#")

    # gene_column_name is supposed to be a column containing
    # primary key of the cancer_susceptibility_genes file

    # check that none of the genes is present multiple times,
    # exit if there is such a gene, report all duplicates
    validate_uniqueness(df, gene_name_column)

    # dictionary of dictionaries to store the input data
    # gene names being the primary keys
    # and for each gene, there is a dictionary
    # with actionability and age
    cancer_susceptibility_genes_dict = {
        row.pop(gene_name_column): row for row in df.to_dicts()
    }

    return (cancer_susceptibility_genes_dict, source)


def generate_predisposition_column_format_description(
    cancer_susceptibility_genes_dict: dict[str, dict[str, str]],
    gene_name_column: str,
    sample_id: str
) -> str:
    """
    Generate text description of columns and values of the cancer susceptibility
    input info. The description will be written in the metadata of the output file.
    """

    # put together column names
    string_format = "[" + gene_name_column + "]"

    unique_inner_keys = sorted(
        {inner_key for inner_dict in cancer_susceptibility_genes_dict.values() for inner_key in inner_dict})

    for key in unique_inner_keys:
        string_format = string_format + "_[" + key + "]"

    description = "# [" + sample_id + \
        "] Column format of gene_predispositions:\n"
    description = description + \
        "# [" + sample_id + "]    a list of character strings of the form:\n"
    description = description + \
        "# [" + sample_id + "]    " + string_format + "\n"
    description = description + "# [" + sample_id + "]    where:\n"

    # for each column of the input data (except for the gene names):
    #  collect possible values,
    #  generate string to print out to description

    # initialize a dictionary where every new key automatically starts with an empty set()
    unique_values_per_key = defaultdict(set)

    # iterate through all inner dictionaries
    for inner_dict in cancer_susceptibility_genes_dict.values():
        for inner_key, inner_value in inner_dict.items():
            # add the value to the set corresponding to this specific inner key
            unique_values_per_key[inner_key].add(inner_value)

    # convert back to a standard dict
    all_keys_values = {k: sorted(list(v))
                       for k, v in unique_values_per_key.items()}

    # add the info to the description string
    for key in sorted(all_keys_values.keys()):
        description = description + \
            "# [" + sample_id + "]        " + key + ": "
        description = description + \
            "[ " + " | ".join(all_keys_values[key]) + " ]" + "\n"

    return description


# TODO: finalise testing from here on


def get_small_variant_calls_file_path(
    sample_id: str,
    config_yaml: str | Path,
    inpred_nomenclature: str | Path,
    root_path: str | Path
) -> Path:
    """
    Return Path to the small variant calls file path.
    """

    small_variants = SmallVariantGenomeVcf(
        config_yaml, inpred_nomenclature, root_path, sample_id)

    return small_variants.path


def get_workflow_id_string(
    config_yaml: str | Path,
    inpred_nomenclature: str | Path,
    root_path: str | Path
) -> str:
    """
    Return workflow id string.
    """
    workflow_output = WorkflowOutput(
        config_yaml, inpred_nomenclature, root_path)

    return workflow_output.workflow_id


def lookup_predisposition_variants(
    cancer_susceptibility_genes_dict: dict[str, dict[str, str]],
    sample_id: str,
    config_yaml: str | Path,
    inpred_nomenclature: str | Path,
    root_path: str | Path
) -> pl.DataFrame:
    """
    Look up predispositions, store all the info from SmallVariantGenomeVcf merged_table and the gene predispositions in the output pl.DataFrame.
    """

    workflow_output = WorkflowOutput(
        config_yaml, inpred_nomenclature, root_path)
    small_variants = SmallVariantGenomeVcf(
        config_yaml, inpred_nomenclature, root_path, sample_id)
    merged_variant_info_pldf = small_variants.merge(workflow_output)

    cancer_susceptibility_genes = list(cancer_susceptibility_genes_dict.keys())

    # select only variants overlapping at least one predisposition gene
    predispositions_basic_pldf = merged_variant_info_pldf.filter(
        pl.col("overlapping_genes").list.set_intersection(
            cancer_susceptibility_genes).list.len() > 0
    )

    # add gene_predisposition column (containing empty list) into the data frame
    predispositions_pldf = predispositions_basic_pldf.with_columns(
        gene_predispositions=pl.lit([]).cast(pl.List(pl.String))
    )

    # transform the data frame into a a dictionary
    predispositions_dict = {
        row.pop("variant_id"): row for row in predispositions_pldf.to_dicts()
    }

    # iterate through all variant_ids
    #   iterate through all overlapping_genes
    #       if overlapping_gene in cancer_susceptibility_genes
    #       then append gene_actionability_age for the overlapping_gene
    for variant_id in predispositions_dict.keys():
        for overlapping_gene in predispositions_dict[variant_id]["overlapping_genes"]:
            if (overlapping_gene in cancer_susceptibility_genes):
                predisposition_str = overlapping_gene
                for k in sorted(cancer_susceptibility_genes):
                    predisposition_str = predisposition_str + "_" + \
                        cancer_susceptibility_genes_dict[overlapping_gene][k]
                predispositions_dict[variant_id]["gene_predisposition"].append(
                    predisposition_str
                )

    # transfrom predispositions_dict to pl.DataFrame

    # transform the nested predispositions_dict into a list of dicts
    # using **fields to unpack the rest of the dictionary values
    rows = [
        {"variant_id": variant_id, **fields}
        for variant_id, fields in predispositions_dict.items()
    ]

    # transform to dataframe
    predispositions_final_pldf = pl.DataFrame(rows)

    return predispositions_final_pldf


def print_predisposition_variants_to_output_file(
    sample_id: str,
    tsoppy_version: str,
    workflow_id: str,
    source: str,
    cancer_susceptibility_genes: Path,
    small_variant_calls: Path,
    length_of_targeted_coding_regions: float,
    tumor_purity: float,
    predisposition_variants: pl.DataFrame,
    predisposition_column_format_description: str,
    output_file: Path,
):
    """
    Print predispositions into the output file.
    """

    header_lines_tsoppy = Template(
        """  # [$sample_id] tsoppy version: $tsoppy_version
"""
    )

    header_lines_workflow = Template(
        """  # [$sample_id] Workflow id: $workflow_id
"""
    )

    header_lines_cancer_susceptibility_data_source = Template(
        """  # [$sample_id] Cancer susceptibility genes are defined in:
# [$sample_id]    file: $cancer_susceptibility_genes
# [$sample_id]    source: $source
"""
    )

    header_lines_small_variant_info_data_source = Template(
        """  # [$sample_id] Small variant calls are defined in: $small_variant_calls
"""
    )

    header_lines_length_of_targeted_coding_regions = Template(
        """  # [$sample_id] Cumulative length of all the targeted coding regions (in millions of bases): $length_of_targeted_coding_regions
"""
    )

    header_lines_tumor_purity = Template(
        """  # [$sample_id] Tumor purity (as a fraction between 0 and 1): $tumor_purity
# [$sample_id]    The tumor purity is provided as an input parameter for tsoppy.
"""
    )

    # open output file for writing
    with open(output_file, "w") as output:
        # print metadata comments to the output file

        # print info about tsoppy and workflow versions
        tsoppy = {"sample_id": sample_id, "tsoppy_version": tsoppy_version}
        output.write(header_lines_tsoppy.safe_substitute(tsoppy))

        workflow = {"sample_id": sample_id, "workflow_id": workflow_id}
        output.write(header_lines_workflow.safe_substitute(workflow))

        # print info about where the susceptibility gene list comes from
        sus_genes_data_source = {
            "sample_id": sample_id,
            "cancer_susceptibility_genes": cancer_susceptibility_genes,
            "source": source,
        }
        output.write(
            header_lines_cancer_susceptibility_data_source.safe_substitute(
                sus_genes_data_source
            )
        )

        # print info about where the small variant info comes from
        small_variant_info = {
            "sample_id": sample_id,
            "small_variant_calls": small_variant_calls,
        }
        output.write(
            header_lines_small_variant_info_data_source.safe_substitute(
                small_variant_info
            )
        )

        targeted_coding_regions = {
            "sample_id": sample_id,
            "length_of_targeted_coding_regions": f"{length_of_targeted_coding_regions:.2f}",
        }
        output.write(
            header_lines_length_of_targeted_coding_regions.safe_substitute(
                targeted_coding_regions
            )
        )

        tumor_purity_info = {
            "sample_id": sample_id,
            "tumor_purity": f"{tumor_purity:.2f}",
        }
        output.write(
            header_lines_tumor_purity.safe_substitute(tumor_purity_info))

        output.write(predisposition_column_format_description)

    # write out the dataframe
    with open(output_file, mode="ab") as output:
        predisposition_variants.write_csv(
            output, include_header=True, separator="\t")

    return


def generate_report(
    output_file: Path,
    tsoppy_version: str,
    sample_id: str,
    config_yaml: str | Path,
    inpred_nomenclature: str | Path,
    root_path: str | Path,
    length_of_targeted_coding_regions: float,
    tumor_purity: float,
    cancer_susceptibility_genes: Path,
    csg_column_list: list[str] = ["Gene", "Actionability", "Age"],
    gene_name_column: str = "Gene",
):
    """
    This function reports variants called by small variant caller that are present
    in the cancer susceptibility genes.
    """

    logger.info("Validate tumor purity range.")
    validate_tumor_purity_range(tumor_purity)

    logger.info(
        f"Load data from the {cancer_susceptibility_genes} input file.")
    cancer_susceptibility_genes_dict, source = (
        load_data_from_cancer_susceptibility_genes_file(
            cancer_susceptibility_genes, csg_column_list, gene_name_column
        )
    )

    logger.info("Generate predisposition column format description.")
    predisposition_column_format_description = generate_predisposition_column_format_description(
        cancer_susceptibility_genes_dict, gene_name_column, sample_id)

    logger.info("Find small variant calls file path.")
    small_variant_calls = get_small_variant_calls_file_path(
        sample_id, config_yaml, inpred_nomenclature, root_path)

    logger.info("Retrieve workflow id string.")
    workflow_id = get_workflow_id_string(
        config_yaml, inpred_nomenclature, root_path)

    logger.info(
        f"Find predisposition_variants in SmallVariantGenomeVcf object of a given sample {sample_id}."
    )
    predisposition_variants = lookup_predisposition_variants(
        cancer_susceptibility_genes_dict, sample_id, config_yaml, inpred_nomenclature, root_path
    )

    logger.info(
        f"Print the predispositions into the {output_file} output file."
    )
    print_predisposition_variants_to_output_file(
        sample_id,
        tsoppy_version,
        workflow_id,
        source,
        cancer_susceptibility_genes,
        small_variant_calls,
        length_of_targeted_coding_regions,
        tumor_purity,
        predisposition_variants,
        predisposition_column_format_description,
        output_file,
    )
