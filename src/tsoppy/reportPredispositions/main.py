"""
This module generates predispositions.
"""


import polars as pl
import logging

logger = logging.getLogger(__name__)

# TODO: lookup_predisposition_variants
# TODO: tests


def load_data_from_cancer_susceptibility_genes_table(cancer_susceptibility_genes):
    """
    Load data from the input table of cancer susceptibility genes.
    """

    # load input data
    df = pl.read_csv(cancer_susceptibility_genes, columns=[
                     "Gene", "Actionability", "Age"], separator="\t")

    # dictionary of dictionaries to store the input data
    # gene names being the primary keys
    # and for each gene, there is a dictionary
    # with actionability and age data
    genes = {
        row.pop("Gene"): row
        for row in df.to_dicts()
    }

    return genes


def lookup_predisposition_variants(genes, small_variant_calls):
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

    return predispositions


def print_header_lines_versions(sample_id, version_string, output_file_handle):
    output_file_handle.write(
        f"# [{sample_id}] Version string: {version_string}\n")
    return


def print_header_lines_cancer_susceptibility_data_source(sample_id, cancer_susceptibility_genes, doi_reference, output_file_handle):
    output_file_handle.write(
        f"# [{sample_id}] Cancer susceptibility genes are defined in:\n")
    output_file_handle.write(
        f"# [{sample_id}] \tfile {cancer_susceptibility_genes}\n")
    output_file_handle.write(
        f"# [{sample_id}] \tarticle {doi_reference}\n")
    return


def print_header_lines_remaining_variant_info_data_source(sample_id, small_variant_calls, output_file_handle):
    output_file_handle.write(
        f"# [{sample_id}] Small variant calls are defined in: {small_variant_calls}\n"
    )
    return


def print_header_lines_length_of_targeted_coding_regions(sample_id, length_of_targeted_coding_regions, output_file_handle):
    # report length of targeted coding regions
    output_file_handle.write(
        f"# [{sample_id}] Cumulative length of all the targeted coding regions (in millions of bases): {length_of_targeted_coding_regions}\n")
    return


def print_header_lines_tumor_purity(sample_id, tumor_purity, output_file_handle):
    # report tumor purity (percentage of tumor cells in the sample)
    output_file_handle.write(
        f"# [{sample_id}] Tumor purity (as a fraction between 0 and 1): {tumor_purity}\n")
    output_file_handle.write(
        f"# [{sample_id}] \tThe tumor purity is provided as an input parameter for tsoppy.\n")
    return


def print_header_lines_gene_predisposition(sample_id, output_file_handle):
    # all the info in the Gene_predisposition column come from the cancer susceptibility genes input file
    output_file_handle.write(
        f"# [{sample_id}] Gene_predisposition column format:\n")
    output_file_handle.write(
        f"# [{sample_id}] \t[Actionability]_[Age]\n")
    output_file_handle.write(
        f"# [{sample_id}] \twhere:\n")
    output_file_handle.write(
        f"# [{sample_id}] \t\tActionability: [ MA-CSG | HA-CSG | SA-CSG ]\n")
    output_file_handle.write(
        f"# [{sample_id}] \t\t\t - MA-CSG = most actionable cancer susceptibility gene, HA-CSG = highly actionable csg, SA-CSG = standardly actionable csg\n")
    output_file_handle.write(
        f"# [{sample_id}] \t\tAge: [ Allages | Age<30 ]\n")
    return


def get_genomic_location(variant_id):
    # variant_id format: chromosome:position:ref>alt
    chromosome, position = variant_id.split(':')[0:2]
    return f"{chromosome}:{position}"


def get_dna_change(variant_id):
    # variant_id format: chromosome:position:ref>alt
    dna_change = variant_id.split(':')[2]
    return dna_change


def print_predisposition_variants_to_output_file(
        sample_id,
        version_string,
        cancer_susceptibility_genes,
        doi_reference,
        small_variant_calls,
        length_of_targeted_coding_regions,
        tumor_purity,
        predisposition_variants,
        output_file):
    """
    Print predispositions into the output file.
    """

    # open output file for writing
    with open(output_file, 'w') as output:

        # print metadata comments to the output file
        print_header_lines_versions(sample_id, version_string, output)
        print_header_lines_cancer_susceptibility_data_source(
            sample_id, cancer_susceptibility_genes, doi_reference, output)
        print_header_lines_remaining_variant_info_data_source(
            sample_id, small_variant_calls, output)
        print_header_lines_length_of_targeted_coding_regions(
            sample_id, length_of_targeted_coding_regions, output)
        print_header_lines_tumor_purity(
            sample_id, tumor_purity, output)
        print_header_lines_gene_predisposition(
            sample_id, output)

    # transform the nested dict into a list of dicts
    # using **fields to unpack the rest of the dictionary values
    rows = [
        {"variant_id": variant_id, **fields}
        for variant_id, fields in predisposition_variants.items()
    ]

    df = pl.DataFrame(rows)

    with open(output_file, mode="ab") as output:
        df.write_csv(output, include_header=True)

    return


def report_predispositions(sample_id, version_string, reference, target_size_coding, tumor_purity, cancer_susceptibility_genes, small_variant_calls, output_file):
    """
    This function reports variants called by small variant caller that are present
    in the cancer susceptibility genes.
    """

    predisposition_variants = dict()
    csg = dict()

    logger.info(f"Load data from the {cancer_susceptibility_genes} table")
    csg = load_data_from_cancer_susceptibility_genes_table(
        cancer_susceptibility_genes)

    logger.info(f"Open the {small_variant_calls} file and iterate through the variants. Store all the variants present in the {cancer_susceptibility_genes} table together with all the info that should be reported into the {predisposition_variants}.")
    predisposition_variants = lookup_predisposition_variants(
        csg, small_variant_calls)

    logger.info(
        f"Print the {predisposition_variants} content into an output file.")
    print_predisposition_variants_to_output_file(
        sample_id, version_string, reference, target_size_coding, tumor_purity, predisposition_variants, output_file)
