
import logging

logger = logging.getLogger()


# TODO: no functionality yet, just documentation
# TODO: should the output file content some patient ID/sample ID to make sure that that one knows patient ID from the output file content? if so, where should the ID be obtained?
# TODO: should output file name be (1) an input parameter or (2) should it be generated automatically, if (2) then where and how
# TODO: will patient age be available at tsoppy run?


def load_data_from_cancer_susceptibility_genes_table(cancer_susceptibility_genes):
    """
    Load data from the input table.
    """
    csg = dict()
    # open input file
    # parse the collumns
    # store the values
    return csg


def lookup_predispositions(csg, small_variant_calls, age):
    """
    Look up predispositions and store all the relevant info.
    """
    predispositions = dict()
    # open small_variant_calls file
    # iterate through all the variants, store variants in the genes present from csg in predispositions, pending on the age
    return predispositions


def print_predispositions_to_output_file(output, predispositions):
    """
    Print predispositions into the output file.
    """
    # open output file for writing
    # iterate through the predispositions, write out


def report_predispositions(cancer_susceptibility_genes, small_variant_calls, age, output):
    """
    This function reports variants called by small variant caller that are present
    in the cancer susceptibility genes. Susceptibility of a gene depends on a patient's age,
    thus the age input parameter.
    """

    predispositions = dict()
    csg = dict()

    logger.info(f"Load data from the {cancer_susceptibility_genes} table")
    csg = load_data_from_cancer_susceptibility_genes_table(
        cancer_susceptibility_genes)

    logger.info(f"Open the {small_variant_calls} file and iterate through the variants. Store all the variants present in the {cancer_susceptibility_genes} table together with all the info that should be reported into the {predispositions}.")
    predispositions = lookup_predispositions(csg, small_variant_calls, age)

    logger.info(f"Print the {predispositions} content into the {output} file.")
    print_predispositions_to_output_file(output, predispositions)
