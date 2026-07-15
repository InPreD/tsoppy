"""
Report predispositions subpackage unit test.
"""

from pathlib import Path

import polars as pl
import pytest

from tsoppy.report_predispositions.main import (
    get_dna_change,
    get_genomic_location,
    load_data_from_cancer_susceptibility_genes_table,
    print_predisposition_variants_to_output_file,
    tumor_purity_range_validation,
    validate_uniqueness,
)

# TODO: test_lookup_predisposition_variants
# TODO: test_generate_report


def test_validate_uniqueness():

    # unique keys
    unique_df = pl.DataFrame(
        {
            "Gene": ["TP53", "BRCA1"],
            "Actionability": ["MA-CSG", "HA-CSG"],
            "Age": ["Age<30", "Allages"],
        }
    )

    assert validate_uniqueness(unique_df, "Gene") is None

    # non_unique keys
    non_unique_df = pl.DataFrame(
        {
            "Gene": ["TP53", "BRCA1", "TP53"],
            "Actionability": ["MA-CSG", "HA-CSG", "SA-CSG"],
            "Age": ["Age<30", "Allages", "Allages"],
        }
    )

    with pytest.raises(ValueError):
        validate_uniqueness(non_unique_df, "Gene")


def test_tumor_purity_range_validation():
    # value in the range
    value = 0.15

    assert tumor_purity_range_validation(value) is None

    # value below 0
    value = -5

    with pytest.raises(ValueError):
        tumor_purity_range_validation(value)

    # value higher than 1
    value = 12

    with pytest.raises(ValueError):
        tumor_purity_range_validation(value)


def test_load_data_from_cancer_susceptibility_genes_table(tmp_path):
    # cases to test:
    # expected input
    # missing column
    # nonexistent file
    # duplicate records for a gene

    column_list = list(("Gene", "Actionability", "Age"))
    gene_name_column = "Gene"

    print(column_list)

    # create temporary directory
    directory = tmp_path / "test_reporPredispositions_main"
    directory.mkdir()

    # write content to the temp file, expected
    test_file_expected = directory / "cancer_susceptibility_genes_expected.csv"
    cancer_susceptibility_genes_table_content_expected = (
        "Gene\tActionability\tAge\nTP53\tHA-CSG\tAge<30\nBRCA1\tMA-CSG\tAllages"
    )
    test_file_expected.write_text(
        cancer_susceptibility_genes_table_content_expected)

    # expected result
    result = dict()
    result["TP53"] = dict()
    result["TP53"]["Actionability"] = "HA-CSG"
    result["TP53"]["Age"] = "Age<30"
    result["BRCA1"] = dict()
    result["BRCA1"]["Actionability"] = "MA-CSG"
    result["BRCA1"]["Age"] = "Allages"

    # test the expected input
    assert (
        load_data_from_cancer_susceptibility_genes_table(
            test_file_expected, column_list, gene_name_column
        )
        == result
    )

    # remove the temporary file test_file_expected
    test_file_expected.unlink()

    # write content to the temp file, missing column
    test_file_missing_column = (
        directory / "cancer_susceptibility_genes_missing_column.csv"
    )
    cancer_susceptibility_genes_table_content_missing_column = (
        "Gene\tActionability\nTP53\tHA-CSG\nBRCA1\tMA-CSG"
    )
    test_file_missing_column.write_text(
        cancer_susceptibility_genes_table_content_missing_column
    )

    # test the input with missing column
    with pytest.raises(pl.exceptions.ColumnNotFoundError):
        load_data_from_cancer_susceptibility_genes_table(
            test_file_missing_column, column_list, gene_name_column
        )

    # remove the temporary file test_file_missing_column
    test_file_missing_column.unlink()

    # write content to the temp file, nonexistent
    test_file_nonexistent = directory / "cancer_susceptibility_genes_nonexistent.csv"

    # test the nonexistent input file
    with pytest.raises(FileNotFoundError):
        load_data_from_cancer_susceptibility_genes_table(
            test_file_nonexistent, column_list, gene_name_column
        )

    # test the input with duplicated genes
    test_file_duplicated = directory / "cancer_susceptibility_genes_duplicated.csv"
    cancer_susceptibility_genes_table_content_duplicated = "Gene\tActionability\tAge\nTP53\tHA-CSG\tAge<30\nBRCA1\tMA-CSG\tAllages\nTP53\tSA-CSG\tAge<30"
    test_file_duplicated.write_text(
        cancer_susceptibility_genes_table_content_duplicated
    )

    with pytest.raises(ValueError):
        load_data_from_cancer_susceptibility_genes_table(
            test_file_duplicated, column_list, gene_name_column
        )

    # remove the temp file test_file_duplicated
    test_file_duplicated.unlink()

    # remove the temp directory
    directory.rmdir()


def test_get_genomic_location():
    variant = "chromosome:position:ref>alt"
    out = "chromosome:position"
    assert get_genomic_location(variant) == out


def test_get_dna_change():
    variant = "chromosome:position:ref>alt"
    out = "ref>alt"
    assert get_dna_change(variant) == out


def test_print_predisposition_variants_to_output_file(tmp_path):

    # create temporary directory
    directory = tmp_path / "test_reportPredispositions_main"
    directory.mkdir()

    sample_id = "IPA0000-X00-Y00-Z00"
    version_string = "VERSION_STRING"

    cancer_susceptibility_genes = directory / "cancer_susceptibility_genes.tsv"
    cancer_susceptibility_genes.write_text(
        "Gene\tActionability\tAge\nTP53\tHA-CSG\tAge<30\nBRCA1\tMA-CSG\tAllages"
    )

    doi_reference = "https://doi.org/10.1016/j.annonc.2022.12.003"

    small_variant_calls = directory / "small_variant_calls.tsv"
    small_variant_calls.write_text("test")

    length_of_targeted_coding_regions = 1.27
    tumor_purity = 0.8

    # create predisposition variants dict
    predisposition_variants = dict()
    predisposition_variants["17:7577144:A>G"] = dict()
    predisposition_variants["17:7577144:A>G"]["sample_id"] = sample_id
    predisposition_variants["17:7577144:A>G"]["gene_symbol"] = "TP53"
    predisposition_variants["17:7577144:A>G"]["ensembl_transcript_id"] = (
        "ENST00000269305.9"
    )
    predisposition_variants["17:7577144:A>G"]["refseq_mrna"] = "NM_000546.6"
    predisposition_variants["17:7577144:A>G"]["genomic_location"] = "17:7577144"
    predisposition_variants["17:7577144:A>G"]["dna_change"] = "A>G"
    predisposition_variants["17:7577144:A>G"]["cdna_change"] = "c.794T>C"
    predisposition_variants["17:7577144:A>G"]["protein_change"] = "p.Leu265Pro"
    predisposition_variants["17:7577144:A>G"]["depth_tumor_dna"] = "648"
    predisposition_variants["17:7577144:A>G"]["af_tumor_dna"] = "0.52"
    predisposition_variants["17:7577144:A>G"]["depth_normal_dna"] = "NA"
    predisposition_variants["17:7577144:A>G"]["af_normal_dna"] = "NA"
    predisposition_variants["17:7577144:A>G"]["depth_tumor_rna"] = "NA"
    predisposition_variants["17:7577144:A>G"]["af_tumor_rna"] = "NA"
    predisposition_variants["17:7577144:A>G"]["tcga_frequency"] = "NA"
    predisposition_variants["17:7577144:A>G"]["icgc_pcawg_occurrence"] = "NA"
    predisposition_variants["17:7577144:A>G"]["gene_predisposition"] = "HA-CSG_Age<30"

    output_file = directory / "predispositions.tsv"

    expected_output_file = directory / "expected_predispositions.tsv"

    metadata = (
        f"# [{sample_id}] Version string: {version_string}\n"
        f"# [{sample_id}] Cancer susceptibility genes are defined in:\n"
        f"# [{sample_id}] \tfile {cancer_susceptibility_genes}\n"
        f"# [{sample_id}] \tarticle {doi_reference}\n"
        f"# [{sample_id}] Small variant calls are defined in: {small_variant_calls}\n"
        f"# [{sample_id}] Cumulative length of all the targeted coding regions (in millions of bases): {length_of_targeted_coding_regions:.2f}\n"
        f"# [{sample_id}] Tumor purity (as a fraction between 0 and 1): {tumor_purity:.2f}\n"
        f"# [{sample_id}] \tThe tumor purity is provided as an input parameter for tsoppy.\n"
        f"# [{sample_id}] Gene_predisposition column format:\n"
        f"# [{sample_id}] \t[Actionability]_[Age]\n"
        f"# [{sample_id}] \twhere:\n"
        f"# [{sample_id}] \t\tActionability: [ MA-CSG | HA-CSG | SA-CSG ]\n"
        f"# [{sample_id}] \t\t\t - MA-CSG = most actionable cancer susceptibility gene, HA-CSG = highly actionable csg, SA-CSG = standardly actionable csg\n"
        f"# [{sample_id}] \t\tAge: [ Allages | Age<30 ]\n"
    )

    header = (
        "sample_id\t"
        "gene_symbol\t"
        "ensembl_transcript_id\t"
        "refseq_mrna\t"
        "genomic_location\t"
        "dna_change\t"
        "cdna_change\t"
        "protein_change\t"
        "depth_tumor_dna\t"
        "af_tumor_dna\t"
        "depth_normal_dna\t"
        "af_normal_dna\t"
        "depth_tumor_rna\t"
        "af_tumor_rna\t"
        "tcga_frequency\t"
        "icgc_pcawg_occurrence\t"
        "gene_predisposition\n"
    )

    variants = (
        f"{sample_id}\t"
        f"TP53\t"
        f"ENST00000269305.9\t"
        f"NM_000546.6\t"
        f"17:7577144\t"
        f"A>G\t"
        f"c.794T>C\t"
        f"p.Leu265Pro\t"
        f"648\t"
        f"0.52\t"
        f"NA\t"
        f"NA\t"
        f"NA\t"
        f"NA\t"
        f"NA\t"
        f"NA\t"
        f"HA-CSG_Age<30\n"
    )

    file_content = metadata + header + variants

    expected_output_file.write_text(file_content)

    print_predisposition_variants_to_output_file(
        sample_id,
        version_string,
        cancer_susceptibility_genes,
        doi_reference,
        small_variant_calls,
        length_of_targeted_coding_regions,
        tumor_purity,
        predisposition_variants,
        output_file,
    )

    assert Path(expected_output_file).read_text() == Path(
        output_file).read_text()
