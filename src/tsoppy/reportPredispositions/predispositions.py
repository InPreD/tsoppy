
import logging
import csv

logger = logging.getLogger(__name__)

# TODO: lookup_predisposition_variants
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


def print_header_lines_copy_number(sample_id, output_file_handle):
    output_file_handle.write(
        f"# [{sample_id}] Gene_CN column format:\n")
    output_file_handle.write(
        f"# [{sample_id}] \t[HC|non-HC]_[Tumor_CN]/[Adjusted_Tumor_CN]_[Normal_CN]\n")
    output_file_handle.write(
        f"# [{sample_id}] \twhere:\n")
    output_file_handle.write(
        f"# [{sample_id}] \t\tCN stands for copy number,\n")
    output_file_handle.write(
        f"# [{sample_id}] \t\tHC stands for high confidence.\n")
    return


def print_header_lines_gene_predisposition(sample_id, output_file_handle):
    output_file_handle.write(
        f"# [{sample_id}] Gene_predisposition column format:\n")
    output_file_handle.write(
        f"# [{sample_id}] \t[category]_[note]_[associated_tumor_type]\n")
    return


def get_genomic_location(variant_id):
    # variant_id format: chromosome:position:ref>alt
    chromosome, position = variant_id.split(':')[0:2]
    return f"{chromosome}:{position}"


def get_dna_change(variant_id):
    # variant_id format: chromosome:position:ref>alt
    dna_change = variant_id.split(':')[2]
    return dna_change


def print_predisposition_variants_to_output_file(sample_id, version_string, cancer_susceptibility_genes, doi_reference, target_size_coding, tumor_purity, predisposition_variants, output_file):
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
            f"# [{sample_id}] Cancer susceptibility genes are defined in:\n")
        output.write(
            f"# [{sample_id}] \tfile {cancer_susceptibility_genes}\n")
        output.write(
            f"# [{sample_id}] \tarticle {doi_reference}\n")
        output.write(
            f"# [{sample_id}] Size of the target coding region (in millions of bases): {target_size_coding}\n")
        output.write(
            f"# [{sample_id}] Specified tumor purity (as a fraction between 0 and 1): {tumor_purity}\n")

        print_header_lines_copy_number(sample_id, output)
        print_header_lines_gene_predisposition(sample_id, output)

        output.write(
            "\t".join(column_names)+"\n")

        # print out body of the output file:
        # iterate through the {predispositions} and for each record write down info for all the columns
        for variant_id in predisposition_variants.keys():
            output_body_line = sample_id + "\t"
            output_body_line += predisposition_variants[variant_id]['gene_symbol'] + "\t"
            output_body_line += predisposition_variants[variant_id]['ensembl_transcript_id'] + "\t"
            output_body_line += predisposition_variants[variant_id]['refseq_mrna'] + "\t"
            output_body_line += get_genomic_location(variant_id) + "\t"
            output_body_line += get_dna_change(variant_id) + "\t"
            output_body_line += predisposition_variants[variant_id]['cdna_change'] + "\t"
            output_body_line += predisposition_variants[variant_id]['protein_change'] + "\t"
            output_body_line += predisposition_variants[variant_id]['depth_tumor_dna'] + "\t"
            output_body_line += predisposition_variants[variant_id]['af_tumor_dna'] + "\t"
            output_body_line += predisposition_variants[variant_id]['depth_normal_dna'] + "\t"
            output_body_line += predisposition_variants[variant_id]['af_normal_dna'] + "\t"
            output_body_line += predisposition_variants[variant_id]['depth_tumor_rna'] + "\t"
            output_body_line += predisposition_variants[variant_id]['af_tumor_rna'] + "\t"
            output_body_line += predisposition_variants[variant_id]['tcga_frequency'] + "\t"
            output_body_line += predisposition_variants[variant_id]['icgc_pcawg_occurrence'] + "\t"
            output_body_line += predisposition_variants[variant_id]['gene_predisposition'] + "\t"
            output_body_line += predisposition_variants[variant_id]['gene_cn'] + "\t"
            output_body_line += predisposition_variants[variant_id]['cpsr_acmg_class'] + "\t"
            output_body_line += predisposition_variants[variant_id]['cpsr_clinvar_class'] + "\t"
            output_body_line += predisposition_variants[variant_id]['cpsr_classification_doc']
            output_body_line += "\n"

            output.write(output_body_line)


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
