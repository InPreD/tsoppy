
import logging
import csv

logger = logging.getLogger()


# TODO: lookup_predisposition_variants
# TODO: print_predisposition_variants_to_output_file
# TODO: tests

def load_data_from_cancer_susceptibility_genes_table(cancer_susceptibility_genes):
    """
    Load data from the input table.
    """
    csg = dict()

    # open the input file
    with open(cancer_susceptibility_genes, mode='r', newline='') as file:
        reader = csv.DictReader(file, delimiter='\t')

        # store all the info in csg dictionary
        for row in reader:
            csg[row['Gene']] = dict()
            csg[row['Gene']]['Actionability'] = row['Actionability']
            csg[row['Gene']]['Age'] = row['Age']
    return csg


def lookup_predisposition_variants(csg, small_variant_calls):
    """
    Look up predispositions and store all the relevant info.
    """
    predispositions = dict()
    # open small_variant_calls file
    # iterate through all the variants, store variants in the genes present from csg in predispositions
    return predispositions


def print_predisposition_variants_to_output_file(sample_id, version_string, reference, target_size_coding, tumor_purity, predisposition_variants, output_file):
    """
    Print predispositions into the output file.
    """

    # define names of columns of the output file
    column_names = [
        "Sample_ID",
        "Gene_symbol",
        "Ensembl_transcript_ID",
        "RefSeq_mRNA",
        "Genomic_location",
        "DNA_change",
        "cDNA_change",
        "Protein_change",
        "Depth_tumor_DNA",
        "AF_tumor_DNA",
        "Depth_normal_DNA",
        "AF_normal_DNA",
        "Depth_tumor_RNA",
        "AF_tumor_RNA",
        "TCGA_frequency",
        "ICGC_PCAWG_occurrence",
        "Gene_predisposition",
        "Gene_CN",
        "CPSR_ACMG_class",
        "CPSR_ClinVar_class",
        "CPSR_classification_doc"]

    # open output file for writing
    with open(output_file, 'w') as output:
        # print out header of the output file
        output.write(
            f"# [{sample_id}] Version string: {version_string}\n")
        output.write(
            f"# [{sample_id}] Variants included in this table are located within one of the cancer predisposition genes listed in {reference}\n")
        output.write(
            f"# [{sample_id}] Size of the target coding region (in millions of bases): {target_size_coding}\n")
        output.write(
            f"# [{sample_id}] Specified tumor purity (as a fraction between 0 and 1): {tumor_purity}\n")
        output.write(
            f"# [{sample_id}] \"Gene_CN\" column format: [High_confidence/non-High_confidence]_[Tumor_CN]/[Adjusted_Tumor_CN]_[Normal_CN]\n")
        output.write(
            "\t".join(column_names)+"\n")

        # print out body of the output file:
        # iterate through the {predispositions} and for each record write down info for all the columns
        for variant_id in predisposition_variants.keys():
            output_body_line = sample_id + "\t"
            output_body_line += predisposition_variants[variant_id]['gene_symbol'] + "\t"
            output_body_line += variant_id
            output_body_line += "\n"
            # TODO: add all the other columns for each predisposition variant

            output.write(f"{output_body_line}")

        # from TSOPPI (user_scripts/libs/05_PCGR_to_variant)interpretation_table.py)
        #
        #   # is the variant located in one of the predisposition genes?
        #   if (SYMBOL in predisposition_gene_values):
        #       with open(
        #           arg_dict["predisposition_output_tsv"], "a") as pot_file:
        #               pot_file.write("\t".join([
        #                   SAMPLE_ID, SYMBOL, ENSEMBL_TRANSCRIPT_ID, REFSEQ_MRNA,
        #                   CHROM + ":" + POS, REF + ">" + ALT, cDNA_change,
        #                   PROTEIN_CHANGE, DP_TUMOR, AF_TUMOR, DP_CONTROL,
        #                   AF_CONTROL, DP_RNA, AF_RNA,
        #                   TCGA_FREQUENCY, ICGC_PCAWG_OCCURRENCE,
        #                   Gene_predisposition, Gene_CN, CPSR_ACMG_class,
        #                   CPSR_ClinVar_class, CPSR_CLASSIFICATION_DOC]) + "\n")


def report_predispositions(cancer_susceptibility_genes, small_variant_calls):
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

    sample_id = "IPA0000-X00-Y00-Z00"
    version_string = "VERSION_STRING"
    reference = "doi_url_esmo_paper_2023"
    target_size_coding = "1.27"
    tumor_purity = "0.8"
    output_file = "/data/test_predisposition_variant_output.csv"
    predisposition_variants["17:7580123:REF>ALT"] = dict()
    predisposition_variants["17:7580123:REF>ALT"]['gene_symbol'] = 'TP53'

    logger.info(
        f"Print the {predisposition_variants} content into an output file.")
    print_predisposition_variants_to_output_file(
        sample_id, version_string, reference, target_size_coding, tumor_purity, predisposition_variants, output_file)
