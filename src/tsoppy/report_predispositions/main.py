"""
This module filters an input list of small variants and outputs a table containing only variants located in cancer-predisposing genes.
"""

import io
import logging
from pathlib import Path

import polars as pl

logger = logging.getLogger(__name__)

# TODO: implement lookup_predisposition_variants


def validate_uniqueness(df: pl.DataFrame, column: str):
    """
    Test that all the values in column 'column' of dataframe 'df' are unique.
    """
    if not df[column].is_unique().all():
        # find the duplicates to make the error message helpful
        duplicates = (
            df[column].filter(df[column].is_duplicated()).unique().to_list()
        )
        # this captures the full traceback automatically
        logger.exception("Data Integrity Validation Failed.")
        raise ValueError(f"Duplicate IDs found in '{column}': {duplicates}")
    else:
        return


def validate_tumor_purity_range(tumor_purity: float):
    """
    Test that tumor purity value is between 0 and 1.
    """
    if not 0 <= tumor_purity <= 1:
        # this captures the full traceback automatically
        logger.exception("Data Integrity Validation Failed.")
        raise ValueError(f"Tumor purity value {tumor_purity} is out of range.")
    else:
        return


def load_data_from_cancer_susceptibility_genes_table(
    file_path: Path, column_list: list[str], gene_name_column: str
) -> dict[str, dict[str, str]]:
    """
    Load data from the input table of cancer susceptibility genes.
    """

    # load input data
    df = pl.read_csv(file_path,
                     columns=column_list, separator="\t")

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

    return cancer_susceptibility_genes_dict


def lookup_predisposition_variants(
    cancer_susceptibility_genes_dict: dict[str, dict[str, str]],
    small_variant_calls: Path,
):
    """
    Look up predispositions and store all the relevant info.
    """
    predispositions = dict()

    #    # define names of columns of the output file
    #    column_names = [
    #        "Sample_ID",
    #        "Gene_symbol",
    #        "Ensembl_transcript_ID",
    #        "RefSeq_mRNA",
    #        "Genomic_location",
    #        "DNA_change",
    #        "cDNA_change",
    #        "Protein_change",
    #        "Depth_tumor_DNA",
    #        "AF_tumor_DNA",
    #        "Depth_normal_DNA",
    #        "AF_normal_DNA",
    #        "Depth_tumor_RNA",
    #        "AF_tumor_RNA",
    #        "TCGA_frequency",
    #        "ICGC_PCAWG_occurrence",
    #        "Gene_predisposition"]

    # open small_variant_calls file
    # iterate through all the variants, store the variants located in the genes in the genes dict in predispositions

    # select only records from the genes in {cancer_susceptibility_genes}

    return predispositions


def print_header_lines_versions(
    sample_id: str, version_string: str, output_file_handle: io.TextIOBase
):
    output_file_handle.write(
        f"# [{sample_id}] Version string: {version_string}\n")
    return


def print_header_lines_cancer_susceptibility_data_source(
    sample_id: str,
    cancer_susceptibility_genes: Path,
    doi_reference: str,
    output_file_handle: io.TextIOBase,
):
    output_file_handle.write(
        f"# [{sample_id}] Cancer susceptibility genes are defined in:\n"
    )
    output_file_handle.write(
        f"# [{sample_id}] \tfile {cancer_susceptibility_genes}\n")
    output_file_handle.write(f"# [{sample_id}] \tarticle {doi_reference}\n")
    return


def print_header_lines_remaining_variant_info_data_source(
    sample_id: str, small_variant_calls: Path, output_file_handle: io.TextIOBase
):
    output_file_handle.write(
        f"# [{sample_id}] Small variant calls are defined in: {small_variant_calls}\n"
    )
    return


def print_header_lines_length_of_targeted_coding_regions(
    sample_id: str,
    length_of_targeted_coding_regions: float,
    output_file_handle: io.TextIOBase,
):
    # report length of targeted coding regions
    output_file_handle.write(
        f"# [{sample_id}] Cumulative length of all the targeted coding regions (in millions of bases): {length_of_targeted_coding_regions:.2f}\n"
    )
    return


def print_header_lines_tumor_purity(
    sample_id: str, tumor_purity: float, output_file_handle: io.TextIOBase
):
    # report tumor purity (percentage of tumor cells in the sample)
    output_file_handle.write(
        f"# [{sample_id}] Tumor purity (as a fraction between 0 and 1): {tumor_purity:.2f}\n"
    )
    output_file_handle.write(
        f"# [{sample_id}] \tThe tumor purity is provided as an input parameter for tsoppy.\n"
    )
    return


def print_header_lines_gene_predisposition(
    sample_id: str, output_file_handle: io.TextIOBase
):
    # all the info in the Gene_predisposition column come from the cancer susceptibility genes input file
    output_file_handle.write(
        f"# [{sample_id}] Gene_predisposition column format:\n")
    output_file_handle.write(f"# [{sample_id}] \t[Actionability]_[Age]\n")
    output_file_handle.write(f"# [{sample_id}] \twhere:\n")
    output_file_handle.write(
        f"# [{sample_id}] \t\tActionability: [ MA-CSG | HA-CSG | SA-CSG ]\n"
    )
    output_file_handle.write(
        f"# [{sample_id}] \t\t\t - MA-CSG = most actionable cancer susceptibility gene, HA-CSG = highly actionable csg, SA-CSG = standardly actionable csg\n"
    )
    output_file_handle.write(
        f"# [{sample_id}] \t\tAge: [ Allages | Age<30 ]\n")
    return


def get_genomic_location(variant_id: str) -> str:
    # variant_id format: chromosome:position:ref>alt
    chromosome, position = variant_id.split(":")[0:2]
    return f"{chromosome}:{position}"


def get_dna_change(variant_id: str) -> str:
    # variant_id format: chromosome:position:ref>alt
    dna_change = variant_id.split(":")[2]
    return dna_change


def print_predisposition_variants_to_output_file(
    sample_id: str,
    version_string: str,
    cancer_susceptibility_genes: Path,
    small_variant_calls: Path,
    length_of_targeted_coding_regions: float,
    tumor_purity: float,
    predisposition_variants: dict[str, dict[str, str | int | float]],
    output_file: Path,
):
    """
    Print predispositions into the output file.
    """

    # open output file for writing
    with open(output_file, "w") as output:
        # print metadata comments to the output file
        print_header_lines_versions(sample_id, version_string, output)
        # -----
        # TODO: get source info from metadata of cancer_susceptibility_genes file and print it out
        # -----
        #  print_header_lines_cancer_susceptibility_data_source(
        #    sample_id, cancer_susceptibility_genes, doi_reference, output
        # )
        print_header_lines_remaining_variant_info_data_source(
            sample_id, small_variant_calls, output
        )
        print_header_lines_length_of_targeted_coding_regions(
            sample_id, length_of_targeted_coding_regions, output
        )
        print_header_lines_tumor_purity(sample_id, tumor_purity, output)
        print_header_lines_gene_predisposition(sample_id, output)

    # transform the nested dict into a list of dicts
    # using **fields to unpack the rest of the dictionary values
    rows = [
        {"variant_id": variant_id, **fields}
        for variant_id, fields in predisposition_variants.items()
    ]

    # transform to dataframe
    df = pl.DataFrame(rows)

    # remove variant_id column from df
    df = df.drop("variant_id")

    # write out the dataframe
    with open(output_file, mode="ab") as output:
        df.write_csv(output, include_header=True, separator="\t")

    return


def generate_report(
    sample_id: str,
    version_string: str,
    length_of_targeted_coding_regions: float,
    tumor_purity: float,
    cancer_susceptibility_genes: Path,
    csg_column_list: list[str],
    gene_name_column: str,
    small_variant_calls: Path,
    output_file: Path,
):
    """
    This function reports variants called by small variant caller that are present
    in the cancer susceptibility genes.
    """

    predisposition_variants = dict()
    cancer_susceptibility_genes_dict = dict()

    tumor_purity_range_validation(tumor_purity)

    logger.info(
        f"Load data from the {cancer_susceptibility_genes} input file.")
    cancer_susceptibility_genes_dict = load_data_from_cancer_susceptibility_genes_table(
        cancer_susceptibility_genes, csg_column_list, gene_name_column
    )

    logger.info(
        f"Open the {small_variant_calls} file and iterate through the variants. Store all the variants present in the {cancer_susceptibility_genes} table together with all the info that should be reported into the {predisposition_variants}."
    )
    predisposition_variants = lookup_predisposition_variants(
        cancer_susceptibility_genes_dict, small_variant_calls
    )

    logger.info(
        f"Print the {predisposition_variants} content into the {output_file} output file."
    )
    print_predisposition_variants_to_output_file(
        sample_id,
        version_string,
        cancer_susceptibility_genes,
        small_variant_calls,
        length_of_targeted_coding_regions,
        tumor_purity,
        predisposition_variants,
        output_file,
    )
