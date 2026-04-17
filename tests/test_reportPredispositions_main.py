"""
Report predispositions subpackage unit test.
"""

import polars as pl
import pytest


from tsoppy.reportPredispositions.main import validate_uniqueness, load_data_from_cancer_susceptibility_genes_table


def test_validate_uniqueness():

    # unique keys
    unique_df = pl.DataFrame(
        {
            "Gene": ["TP53", "BRCA1"],
            "Actionability": ["MA-CSG", "HA-CSG"],
            "Age": ["Age<30", "Allages"]
        }
    )

    assert validate_uniqueness(unique_df, "Gene") is None

    # non_unique keys

    non_unique_df = pl.DataFrame(
        {
            "Gene": ["TP53", "BRCA1", "TP53"],
            "Actionability": ["MA-CSG", "HA-CSG", "SA-CSG"],
            "Age": ["Age<30", "Allages", "Allages"]
        }
    )

    with pytest.raises(ValueError):
        validate_uniqueness(non_unique_df, "Gene")


def test_load_data_from_cancer_susceptibility_genes_table(tmp_path):
    # cases to test:
    # expected input
    # missing column
    # nonexistent file
    # duplicate records for a gene

    column_list = list(("Gene", "Actionability", "Age"))

    print(column_list)

    # create temporary directory
    directory = tmp_path / 'test_reporPredispositions_main'
    directory.mkdir()

    # write content to the temp file, expected
    test_file_expected = directory / 'cancer_susceptibility_genes_expected.csv'
    cancer_susceptibility_genes_table_content_expected = "Gene\tActionability\tAge\nTP53\tHA-CSG\tAge<30\nBRCA1\tMA-CSG\tAllages"
    test_file_expected.write_text(
        cancer_susceptibility_genes_table_content_expected)

    # expected result
    result = dict()
    result['TP53'] = dict()
    result['TP53']['Actionability'] = 'HA-CSG'
    result['TP53']['Age'] = 'Age<30'
    result['BRCA1'] = dict()
    result['BRCA1']['Actionability'] = 'MA-CSG'
    result['BRCA1']['Age'] = 'Allages'

    # test the expected input
    assert load_data_from_cancer_susceptibility_genes_table(
        test_file_expected, column_list) == result

    # write content to the temp file, missing column
    test_file_missing_column = directory / \
        'cancer_susceptibility_genes_missing_column.csv'
    cancer_susceptibility_genes_table_content_missing_column = "Gene\tActionability\nTP53\tHA-CSG\nBRCA1\tMA-CSG"
    test_file_missing_column.write_text(
        cancer_susceptibility_genes_table_content_missing_column)

    # test the input with missing column
    with pytest.raises(pl.exceptions.ColumnNotFoundError):
        load_data_from_cancer_susceptibility_genes_table(
            test_file_missing_column, column_list)

    # write content to the temp file, nonexistent
    test_file_nonexistent = directory / 'cancer_susceptibility_genes_nonexistent.csv'

    # test the nonexistent input file
    with pytest.raises(FileNotFoundError):
        load_data_from_cancer_susceptibility_genes_table(
            test_file_nonexistent, column_list)

    # test the input with duplicated genes
    test_file_duplicated = directory / 'cancer_susceptibility_genes_duplicated.csv'
    cancer_susceptibility_genes_table_content_duplicated = "Gene\tActionability\tAge\nTP53\tHA-CSG\tAge<30\nBRCA1\tMA-CSG\tAllages\nTP53\tSA-CSG\tAge<30"
    test_file_duplicated.write_text(
        cancer_susceptibility_genes_table_content_duplicated)

    with pytest.raises(ValueError):
        load_data_from_cancer_susceptibility_genes_table(
            test_file_duplicated, column_list)


"""
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

Unit test for the report predisposition function of the subpackage.

    assert True
"""
