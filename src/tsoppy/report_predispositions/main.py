"""
This module filters an input list of small variants and outputs a table containing only the variants located in cancer-predisposing genes.
"""

import logging
import re
from collections import defaultdict
from pathlib import Path
from polars.exceptions import ColumnNotFoundError
from mako.template import Template as makoTemplate

import polars as pl
import polars.selectors as cs
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

    logger.debug(
        f"Check that column \"{column}\" is a primary key of the data frame.")

    if not df[column].is_unique().all():
        # find the duplicates to make the error message helpful
        duplicates = df[column].filter(
            df[column].is_duplicated()).unique().to_list()
        # this captures the full traceback automatically
        logger.exception(
            f"Data Integrity Validation Failed. Duplicate IDs found in column '{column}': {duplicates}.")
        raise ValueError()
    else:
        return


def validate_tumor_purity_range(tumor_purity: float):
    """
    Test that tumor purity value is between 0 and 1.
    """

    logger.debug(
        f"Test that tumor purity value {tumor_purity} is between 0 and 1.")

    if not 0 <= tumor_purity <= 1:
        # this captures the full traceback automatically
        logger.exception(
            f"Data Integrity Validation Failed. Tumor purity value {tumor_purity} is out of range.")
        raise ValueError()
    else:
        return


def check_gene_name_column_presence_in_column_list(gene_name_column: str, column_list: list[str]):
    """
    Check that gene_name_column is present in the provided column_list.
    """

    logger.debug(
        f"Check that column \"{gene_name_column}\" is present in column list {column_list}.")

    if not gene_name_column in column_list:
        logger.exception(
            f"Data Integrity Validation Failed. Column {gene_name_column} is not present in provided column list {column_list}.")
        raise ValueError()
    else:
        return


def parse_metadata_to_get_source(file_path: Path
                                 ) -> str:
    """
    Parse metadata of the file to get source of the information stored in the file.
    """

    logger.debug(
        f"Parsing metadata of {file_path} to get source of the information stored in the file.")

    # try to get metadata about source of the info
    source = "NA"

    with open(file_path, "r") as file:
        for line in file:
            if re.match(r"# source:", line):
                source = line.replace("# source:", "", 1)
                source = source.strip()
                return source

    logger.exception(
        f"Data Integrity Validation Failed. Source information not found in metadata of {file_path}.")
    raise ValueError()


def check_columns_presence_in_file(column_list: list[str], file_path: Path):
    """
    Check that all requested columns are present in the input file.
    """

    logger.debug(
        f"Check that all requested columns from the column_list {column_list} are present in the file {file_path}.")

    lf = pl.scan_csv(file_path, separator="\t", comment_prefix="#")

    available_columns = lf.collect_schema().names()
    missing_columns = [
        col for col in column_list if col not in available_columns]

    if missing_columns:
        logger.error(
            f"Failed to read the file {file_path} because a required column is missing, missing columns: {missing_columns}.")
        raise ColumnNotFoundError()
    else:
        return


def load_data_from_cancer_susceptibility_genes_file(
    file_path: Path, column_list: list[str], gene_name_column: str
) -> tuple[dict[str, dict[str, str]], str]:
    """
    Load data from the input cancer susceptibility genes file.
    """

    logger.debug(
        f"Load data from the cancer susceptibility genes file {file_path}.")

    # check that gene_name_column is present in the column_list
    check_gene_name_column_presence_in_column_list(
        gene_name_column, column_list)

    # get intormation about the source of data presented in the input file
    source = parse_metadata_to_get_source(file_path)

    # check that all columns from the column list are present in the file
    check_columns_presence_in_file(column_list, file_path)

    logger.debug(f"Read tabular content of the file {file_path}.")
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
    config_yaml: Path,
    inpred_nomenclature: Path,
    root_path: Path
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
    config_yaml: Path,
    inpred_nomenclature: Path,
    root_path: Path
) -> pl.DataFrame:
    """
    Look up predispositions, store all the info from SmallVariantGenomeVcf merged_table and the gene predispositions in the output pl.DataFrame.
    """

    logger.debug(
        "Look up overlaps between the small variant calls and the cancer susceptibility genes.")

    logger.debug("Read in merged_table of small variants.")
    workflow_output = WorkflowOutput(
        config_yaml, inpred_nomenclature, root_path)
    small_variants = SmallVariantGenomeVcf(
        config_yaml, inpred_nomenclature, root_path, sample_id)
    merged_variant_info_df = small_variants.merge(workflow_output)

    cancer_susceptibility_genes = list(cancer_susceptibility_genes_dict.keys())

    logger.debug(
        "Select only those small variants which overlap with at least one predispositon gene.")
    # select only variants overlapping at least one predisposition gene
    predispositions_basic_df = merged_variant_info_df.filter(
        pl.col("overlapping_genes").list.set_intersection(
            cancer_susceptibility_genes).list.len() > 0
    )

    if predispositions_basic_df.is_empty():
        logger.debug(
            "There are no small variants overlapping a predisposition gene.")
        return pl.DataFrame(
            {
                "variant_ID": [],
                "variant_type": [],
                "Illumina_variant_class": [],
                "overlapping_genes": [],
                "AD": [],
                "DP": [],
                "AF": [],
                "gene_predispositions": [],
            }
        )

    logger.debug(
        "Fill predisposition gene info for every small variant overlapping the gene.")

    # add gene_predisposition column (containing empty list) into the data frame
    predispositions_df = predispositions_basic_df.with_columns(
        gene_predispositions=pl.lit([]).cast(pl.List(pl.String))
    )

    # transform the data frame into a a dictionary
    predispositions_dict = {
        row.pop("variant_ID"): row for row in predispositions_df.to_dicts()
    }

    # iterate through all variant_ids
    #   iterate through all overlapping_genes
    #       if overlapping_gene in cancer_susceptibility_genes
    #       then append _ separated list of all values for the overlapping gene
    #       _geneactionability_age as of now
    for variant_ID in predispositions_dict.keys():
        for overlapping_gene in predispositions_dict[variant_ID]["overlapping_genes"]:
            if (overlapping_gene in cancer_susceptibility_genes):
                predisposition_str = overlapping_gene

                for k in sorted(cancer_susceptibility_genes_dict[overlapping_gene]):
                    predisposition_str = predisposition_str + "_" + \
                        cancer_susceptibility_genes_dict[overlapping_gene][k]
                predispositions_dict[variant_ID]["gene_predispositions"].append(
                    predisposition_str
                )

    # transfrom predispositions_dict to pl.DataFrame

    # transform the nested predispositions_dict into a list of dicts
    # using **fields to unpack the rest of the dictionary values
    rows = [
        {"variant_ID": variant_ID, **fields}
        for variant_ID, fields in predispositions_dict.items()
    ]

    # transform to dataframe
    predispositions_final_df = pl.DataFrame(rows)

    return predispositions_final_df


def format_predisposition_variants_before_printing(
    sample_id: str,
    predisposition_variants: pl.DataFrame,
) -> pl.DataFrame:
    """
    Format predisposition_variants data frame for printing.
    """

    logger.debug("Formatting predisposition_variants data frame for printing.")

    # if the input dataframe is empty,
    # return empty dataframe of appropriate structure
    if predisposition_variants.is_empty():
        return pl.DataFrame(
            {
                "sample_id": [],
                "overlapping_genes": [],
                "genomic_location": [],
                "DNA_change": [],
                "depth_tumor_DNA": [],
                "AF_tumor_DNA": [],
                "gene_predispositions": [],
            }
        )

    # continue from here on with non-empty dataframe

    # add sample_id clolumn
    predisposition_variants = predisposition_variants.with_columns(
        pl.lit(sample_id).alias("sample_id")
    )

    # split variant_id to genomic_location and DNA_change
    # - split variant_id to chromosome, position and DNA_change
    predisposition_variants = predisposition_variants.with_columns(
        pl.col("variant_ID")
        .str.split_exact(":", 2)
        .struct.rename_fields(["chromosome", "position", "DNA_change"])
        .alias("split_data")
    ).unnest("split_data")
    # - join chromosome and position to genomic_location, drop the partial columns
    predisposition_variants = predisposition_variants.with_columns(
        pl.concat_str(["chromosome", "position"], separator=":")
        .alias("genomic_location")
    ).drop("chromosome", "position")

    # format DP to depth_tumor_DNA: simple renamed copy
    predisposition_variants = predisposition_variants.with_columns(
        pl.col("DP").alias("depth_tumor_DNA")
    )

    # format AF to AF_tumor_DNA: select only the first element of the array
    # round to 4 decimal places
    predisposition_variants = predisposition_variants.with_columns(
        pl.concat_list("AF").list.get(0)
        .alias("AF_tumor_DNA")
    )

    # define output column order
    col_order = [
        "sample_id",
        "overlapping_genes",
        "genomic_location",
        "DNA_change",
        "depth_tumor_DNA",
        "AF_tumor_DNA",
        "gene_predispositions"
    ]

    # make the dataframe printable to a csv file

    # join values of complex columns into comma separated strings
    # columns affected: overlapping_genes, AD, AF, gene_predispositions
    # format ints to strings
    # columns affected: depth_tumor_DNA
    # format floats to 4decimal precision strings / format floats to strings
    # columns affected: AF_tumor_DNA
    predisposition_variants_str = predisposition_variants.with_columns(
        cs.by_dtype(pl.List).list.eval(
            pl.element().cast(pl.String)).list.join(", "),
        cs.integer().cast(pl.String),
        cs.float().map_elements(repr, return_dtype=pl.String),
        # cs.float().cast(pl.String),
        # cs.float().map_elements(
        #    lambda x: f"{x:.4f}" if x is not None else None,
        #    return_dtype=pl.String
    ).select(col_order)

    return predisposition_variants_str


def print_output_file(
    sample_id: str,
    tsoppy_version: str,
    workflow_id: str,
    source: str,
    cancer_susceptibility_genes: Path,
    small_variant_calls: Path,
    coding_region_size_in_megabases: float,
    tumor_purity: float,
    cancer_susceptibility_genes_dict: dict[str, dict[str, str]],
    gene_name_column: str,
    predisposition_variants: pl.DataFrame,
    output_file: Path,
):
    """
    Generate output file from a mako template.
    """

    logger.debug("Generate output file from a mako template.")

    # TODO: add columns depth_normal_DNA, AF_normal_DNA, depth_tumor_RNA, AF_tumor_RNA, assign appropriate values
    # TODO: this should be happening in the lookup function maybe?

    my_template = makoTemplate(
        filename='/workspaces/tsoppy/src/tsoppy/report_predispositions/templates/predispositions.mako')

    # generate list of columns present in the cancer susceptibility genes file
    column_list = sorted(
        {inner_key for inner_dict in cancer_susceptibility_genes_dict.values() for inner_key in inner_dict})

    # insert gene_name_column to the beginning of the list
    column_list.insert(0, gene_name_column)

    # for each column of the cancer susceptibility genes data (except for the gene names):
    #  collect possible values (unique)

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

    predisposition_variants_str = format_predisposition_variants_before_printing(
        sample_id, predisposition_variants)

    rendered_content = my_template.render(
        sample_id=sample_id,
        tsoppy_version=tsoppy_version,
        workflow_id=workflow_id,
        cancer_susceptibility_genes=cancer_susceptibility_genes,
        source=source,
        small_variant_calls=small_variant_calls,
        coding_region_size_in_megabases=coding_region_size_in_megabases,
        tumor_purity=tumor_purity,
        column_list=column_list,
        all_columns_values=all_keys_values,
        predispositions=predisposition_variants_str,
    )

    with open(output_file, mode="w", encoding="utf-8") as output:
        output.write(rendered_content)

    return


def generate_report(
    output_file: Path,
    tsoppy_version: str,
    sample_id: str,
    config_yaml: Path,
    inpred_nomenclature: Path,
    root_path: Path,
    tumor_purity: float,
    cancer_susceptibility_genes: Path,
    csg_column_list: list[str],
    gene_name_column: str,
):
    """
    This function reports variants called by small variant caller that are present
    in the cancer susceptibility genes.
    """

    logger.info(
        "Report small variants presents in the cancer susceptibility genes.")

    logger.debug("Validate tumor purity range.")
    validate_tumor_purity_range(tumor_purity)

    logger.debug(
        f"Load data from the {cancer_susceptibility_genes} input file.")
    cancer_susceptibility_genes_dict, source = (
        load_data_from_cancer_susceptibility_genes_file(
            cancer_susceptibility_genes, csg_column_list, gene_name_column
        )
    )

    # TODO: use functionality from the other module
    # "Cummulative length of all coding regions targeted in the TSO500 experiment (in millions of bases)."
    logger.debug(
        "Retrieve \"Coding Region Size in Megabases\" from the CombinedVariant Output file.")
    coding_region_size_in_megabases = 1.11

    logger.debug("Find small variant calls file path.")
    small_variant_calls = get_small_variant_calls_file_path(
        sample_id, config_yaml, inpred_nomenclature, root_path)

    logger.debug("Retrieve workflow id string.")
    workflow_id = get_workflow_id_string(
        config_yaml, inpred_nomenclature, root_path)

    logger.debug(
        f"Find predisposition_variants in SmallVariantGenomeVcf object of a given sample {sample_id}."
    )
    predisposition_variants = lookup_predisposition_variants(
        cancer_susceptibility_genes_dict, sample_id, config_yaml, inpred_nomenclature, root_path
    )

    logger.debug(
        f"Print the predispositions into the {output_file} output file."
    )
    print_output_file(
        sample_id,
        tsoppy_version,
        workflow_id,
        source,
        cancer_susceptibility_genes,
        small_variant_calls,
        coding_region_size_in_megabases,
        tumor_purity,
        cancer_susceptibility_genes_dict,
        gene_name_column,
        predisposition_variants,
        output_file,
    )
