"""
Report predispositions subpackage unit test.
"""

import os

from tsoppy.reportPredispositions.main import load_data_from_cancer_susceptibility_genes_table, report_predispositions


def test_load_data_from_cancer_susceptibility_genes_table(tmp_path):

    # create temporary file
    directory = tmp_path / 'tsoppy_test'
    directory.mkdir()
    test_file = directory / 'cancer_susceptibility_genes.csv'

    # write content to the temp file
    cancer_susceptibility_genes_table_content = "Gene\tActionability\tAge\nTP53\tHA-CSG\tAge<30\n"
    test_file.write_text(cancer_susceptibility_genes_table_content)

    # create result
    result = dict()
    result['TP53'] = dict()
    result['TP53']['Actionability'] = 'HA-CSG'
    result['TP53']['Age'] = 'Age<30'

    # pass the file to the function
    assert load_data_from_cancer_susceptibility_genes_table(
        test_file) == result


def test_print_predisposition_variants_to_output_file():

    # create predisposition variants dict
    predisposition_variants = dict()
    predisposition_variants['17:7577144:A>G'] = dict()
    predisposition_variants['17:7577144:A>G']['gene_symbol'] = 'TP53'
    predisposition_variants['17:7577144:A>G']['ensembl_transcript_id'] = 'ENST00000269305.9'
    predisposition_variants['17:7577144:A>G']['refseq_mrna'] = 'NM_000546.6'
    predisposition_variants['17:7577144:A>G']['cdna_change'] = 'c.794T>C'
    predisposition_variants['17:7577144:A>G']['protein_change'] = 'p.Leu265Pro'
    predisposition_variants['17:7577144:A>G']['depth_tumor_dna'] = '648'
    predisposition_variants['17:7577144:A>G']['af_tumor_dna'] = '0.52'
    predisposition_variants['17:7577144:A>G']['depth_normal_dna'] = 'NA'
    predisposition_variants['17:7577144:A>G']['af_normal_dna'] = 'NA'
    predisposition_variants['17:7577144:A>G']['depth_tumor_rna'] = 'NA'
    predisposition_variants['17:7577144:A>G']['af_tumor_rna'] = 'NA'
    predisposition_variants['17:7577144:A>G']['tcga_frequency'] = 'NA'
    predisposition_variants['17:7577144:A>G']['icgc_pcawg_occurrence'] = 'NA'
    predisposition_variants['17:7577144:A>G']['gene_predisposition'] = '4_c_NS'
    predisposition_variants['17:7577144:A>G']['gene_cn'] = '[non-HC]_1.94/1.98_NA'
    predisposition_variants['17:7577144:A>G']['cpsr_acmg_class'] = 'Likely_Pathogenic'
    predisposition_variants['17:7577144:A>G']['cpsr_clinvar_class'] = 'Pathogenic'
    predisposition_variants['17:7577144:A>G']['cpsr_classification_doc'] = 'description'

    assert True

# TODO
#   sample_id = "IPA0000-X00-Y00-Z00"
#    version_string = "VERSION_STRING"
#    reference = "doi_url_esmo_paper_2023"
#    target_size_coding = "1.27"
#    tumor_purity = "0.8"
#    output_file = "/data/test_predisposition_variant_output.csv"


def test_report_predispositions():
    """
    Unit test for the report predisposition function of the subpackage.
    """
    assert True
