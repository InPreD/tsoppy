"""
This module filters an input list of small variants and outputs a table containing only variants located in cancer-predisposing genes.
"""

import logging
import re
from pathlib import Path
from string import Template

import polars as pl

logger = logging.getLogger(__name__)

# TODO: implement lookup_predisposition_variants


def validate_uniqueness(df: pl.DataFrame, column: str):
    """
    Test that all the values in column 'column' of dataframe 'df' are unique.
    """
    if not df[column].is_unique().all():
        # find the duplicates to make the error message helpful
        duplicates = df[column].filter(df[column].is_duplicated()).unique().to_list()
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
) -> tuple[dict[str, dict[str, str]], str]:
    """
    Load data from the input table of cancer susceptibility genes.
    """

    # load metadata about source of the info
    with open(file_path, "r") as file:
        for line in file:
            if re.match(r"# source:", line):
                source = line.replace("# source:", "", 1)
                source = source.strip()

    # load input data
    df = pl.read_csv(file_path, columns=column_list, separator="\t")

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

    return cancer_susceptibility_genes_dict, source


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
    source: str,
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

    header_lines_versions = Template(
        """# [$sample_id] Version string: $version_string
"""
    )

    header_lines_cancer_susceptibility_data_source = Template(
        """# [$sample_id] Cancer susceptibility genes are defined in:
# [$sample_id]    file: $cancer_susceptibility_genes
# [$sample_id]    source: $source
"""
    )

    header_lines_small_variant_info_data_source = Template(
        """# [$sample_id] Small variant calls are defined in: $small_variant_calls
"""
    )

    header_lines_length_of_targeted_coding_regions = Template(
        """# [$sample_id] Cumulative length of all the targeted coding regions (in millions of bases): $length_of_targeted_coding_regions}
"""
    )

    header_lines_tumor_purity = Template(
        """# [$sample_id] Tumor purity (as a fraction between 0 and 1): $tumor_purity
# [$sample_id]    The tumor purity is provided as an input parameter for tsoppy.
"""
    )

    header_lines_gene_predisposition = Template(
        """# [$sample_id] Gene_predisposition column format:
# [$sample_id]    [Actionability]_[Age]
# [$sample_id]    where:
# [$sample_id]        Actionability: [ MA-CSG | HA-CSG | SA-CSG ]
# [$sample_id]            - MA-CSG = most actionable cancer susceptibility gene, HA-CSG = highly actionable csg, SA-CSG = standardly actionable csg
# [$sample_id]        Age: [ Allages | Age<30 ]
"""
    )
    # all the info in the Gene_predisposition column come from the cancer susceptibility genes input file

    # open output file for writing
    with open(output_file, "w") as output:
        # print metadata comments to the output file

        # print info about tool versions
        versions = {"sample_id": sample_id, "version_string": version_string}
        output.write(header_lines_versions.safe_substitute(versions))

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
        output.write(header_lines_tumor_purity.safe_substitute(tumor_purity_info))

        gene_predisposition_info = {"sample_id": sample_id}
        header_lines_gene_predisposition.safe_substitute(gene_predisposition_info)

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

    tumor_purity_range_validation(tumor_purity)

    logger.info(f"Load data from the {cancer_susceptibility_genes} input file.")
    cancer_susceptibility_genes_dict, source = (
        load_data_from_cancer_susceptibility_genes_table(
            cancer_susceptibility_genes, csg_column_list, gene_name_column
        )
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
        source,
        cancer_susceptibility_genes,
        small_variant_calls,
        length_of_targeted_coding_regions,
        tumor_purity,
        predisposition_variants,
        output_file,
    )
