
import logging

logger = logging.getLogger()


def report_predispositions(cancer_susceptibility_genes, small_variant_calls, age, output):
    """
    This function reports variants called by small variant caller that are present 
    in the cancer susceptibility genes. Susceptibility of a gene depends on a patients age group, 
    thus age input parameter.
    """

    predispositions = dict()

    logger.info(f"Load data from the {cancer_susceptibility_genes} table")
    logger.info(f"Open the {small_variant_calls} file and iterate through the variants. Store all the variants present in the {cancer_susceptibility_genes} table together with all the info that should be reported into the {predispositions}.")
    logger.info(f"Print the {predispositions} content into the {output} file.")
