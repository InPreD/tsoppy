"""
Report predispositions subpackage unit test.
"""

from contextlib import nullcontext
from pathlib import Path
from polars.exceptions import (ColumnNotFoundError, NoDataError)

import polars as pl
import pytest

from tsoppy.report_predispositions.main import (
    check_columns_presence_in_file,
    check_file_existence_and_content_presence,
    check_gene_name_column_presence_in_column_list,
    generate_predisposition_column_format_description,
    generate_report,
    get_dna_change,
    get_genomic_location,
    get_small_variant_calls_file_path,
    load_data_from_cancer_susceptibility_genes_file,
    lookup_predisposition_variants,
    parse_metadata_to_get_source,
    print_predisposition_variants_to_output_file,
    validate_tumor_purity_range,
    validate_uniqueness,
)


@pytest.mark.parametrize(
    "inputs, exception, want",
    [
        (
            (
                pl.DataFrame(
                    {
                        "Gene": [
                            "TP53",
                            "BRCA1"
                        ],
                        "Actionability": [
                            "MA-CSG",
                            "HA-CSG"
                        ],
                        "Age": [
                            "Age<30",
                            "Allages"
                        ],
                    }
                ),
                "Gene",
            ),
            nullcontext(),
            None,
        ),
        (
            (
                pl.DataFrame(
                    {
                        "Gene": [
                            "TP53",
                            "BRCA1",
                            "TP53"
                        ],
                        "Actionability": [
                            "MA-CSG",
                            "HA-CSG",
                            "SA-CSG"
                        ],
                        "Age": [
                            "Age<30",
                            "Allages",
                            "Allages"
                        ],
                    }
                ),
                "Gene",
            ),
            pytest.raises(ValueError),
            None,
        )
    ],
)
def test_validate_uniqueness(inputs, exception, want):
    with exception:
        got = validate_uniqueness(*inputs)
        assert got == want


@pytest.mark.parametrize(
    "inputs, exception, want",
    [
        (
            (
                0.15,
            ),
            nullcontext(),
            None,
        ),
        (
            (
                -5,
            ),
            pytest.raises(ValueError),
            None,
        ),
        (
            (
                12,
            ),
            pytest.raises(ValueError),
            None
        ),
    ],
)
def test_validate_tumor_purity_range(inputs, exception, want):
    with exception:
        got = validate_tumor_purity_range(*inputs)
        assert got == want


@pytest.mark.parametrize(
    "inputs, exception, want",
    [
        (
            (
                "Gene",
                [
                    "Gene",
                    "Actionability",
                    "Age"
                ]
            ),
            nullcontext(),
            None,
        ),
        (
            (
                "Gene",
                [
                    "Actionability",
                    "Age"
                ]
            ),
            pytest.raises(ValueError),
            None,
        ),
    ],
)
def test_check_gene_name_column_presence_in_column_list(inputs, exception, want):
    with exception:
        got = check_gene_name_column_presence_in_column_list(*inputs)
        assert got == want


@pytest.mark.parametrize(
    "inputs, exception, want",
    # inputs: filename, content
    [
        (
            (
                "nonexistent.csv",
                None
            ),
            pytest.raises(FileNotFoundError),
            None
        ),
        (
            (
                "empty.csv",
                ""
            ),
            pytest.raises(NoDataError),
            None
        ),
        (
            (
                "nonempty.csv",
                "content"
            ),
            nullcontext(),
            None
        ),
    ]
)
def test_check_file_existence_and_content_presence(tmp_path, inputs, exception, want):
    directory = tmp_path / "test_reportPredispositions_main" / "test_check_file"
    directory.mkdir(parents=True, exist_ok=True)

    if not inputs[0] == "nonexistent.csv":
        test_file = directory / inputs[0]
        test_file.write_text(inputs[1], encoding="utf-8")
    else:
        test_file = directory / inputs[0]

    with exception:
        got = check_file_existence_and_content_presence(test_file)
        assert got == want

    if not inputs[0] == "nonexistent.csv":
        test_file.unlink()

    # remove the temp directory
    directory.rmdir()


@pytest.mark.parametrize(
    "inputs, exception, want",
    # inputs: (file_name, content)
    [
        (
            ("contains_source_info.csv", "# source: source_url"),
            nullcontext(),
            "source_url"
        ),
        (
            ("no_source_info.csv", "# irrelevant info"),
            pytest.raises(ValueError),
            None
        ),
    ]
)
def test_parse_metadata_to_get_source(tmp_path, inputs, exception, want):

    directory = tmp_path / "test_reportPredispositions_main" / "test_parse_metadata"
    directory.mkdir(parents=True, exist_ok=True)

    test_file = directory / inputs[0]
    test_file.write_text(inputs[1], encoding="utf-8")

    with exception:
        got = parse_metadata_to_get_source(test_file)
        assert got == want

    test_file.unlink()

    # remove the temp directory
    directory.rmdir()


@pytest.mark.parametrize(
    "inputs, exception, want",
    # inputs: filename, content, column_list
    [
        (
            (
                "missing_column.csv",
                "# source: doi_url\nGene\tActionability\tAge\nTP53\tHA-CSG\tAge<30\nBRCA1\tMA-CSG\tAllages",
                [
                    "Gene",
                    "Actionability",
                    "Age",
                    "X"
                ],
            ),
            pytest.raises(ColumnNotFoundError),
            None
        ),
        (
            (
                "expected.csv",
                "# source: doi_url\nGene\tActionability\tAge\nTP53\tHA-CSG\tAge<30\nBRCA1\tMA-CSG\tAllages",
                [
                    "Gene",
                    "Actionability",
                    "Age"
                ],
            ),
            nullcontext(),
            None
        ),
    ]
)
def test_check_columns_presence_in_file(tmp_path, inputs, exception, want):

    # create temporary directory
    directory = tmp_path / "test_reporPredispositions_main" / "test_check_columns"
    directory.mkdir(parents=True, exist_ok=True)

    # write content to the temp file
    test_file = directory / inputs[0]
    test_file.write_text(inputs[1], encoding="utf-8")

    with exception:
        got = check_columns_presence_in_file(
            inputs[2], test_file)
        assert got == want

    # remove the temporary file test_file
    test_file.unlink()

    # remove the temp directory
    directory.rmdir()


@pytest.mark.parametrize(
    "inputs, want",
    # inputs: filename, content, column_list, gene_name_column
    [
        (
            (
                "cancer_susceptibility_genes_expected.csv",
                "# source: doi_url\nGene\tActionability\tAge\nTP53\tHA-CSG\tAge<30\nBRCA1\tMA-CSG\tAllages",
                [
                    "Gene",
                    "Actionability",
                    "Age"
                ],
                "Gene"
            ),
            (
                {
                    "TP53": {
                        "Actionability": "HA-CSG",
                        "Age": "Age<30",
                    },
                    "BRCA1": {
                        "Actionability": "MA-CSG",
                        "Age": "Allages",
                    },
                },
                "doi_url"
            ),
        ),
    ]
)
def test_load_data_from_cancer_susceptibility_genes_file(tmp_path, inputs, want):

    # create temporary directory
    directory = tmp_path / "test_reporPredispositions_main" / "test_load_data"
    directory.mkdir(parents=True, exist_ok=True)

    # write content to the temp file
    test_file = directory / inputs[0]
    test_file.write_text(inputs[1], encoding="utf-8")

    got = load_data_from_cancer_susceptibility_genes_file(
        test_file, *inputs[2:4])
    assert got == want

    # remove the temporary file test_file
    test_file.unlink()

    # remove the temp directory
    directory.rmdir()


@pytest.mark.parametrize(
    "inputs, want",
    [
        (
            (
                {
                    "TP53": {
                        "Actionability": "HA-CSG",
                        "Age": "Age<30",
                    },
                    "BRCA1": {
                        "Actionability": "MA-CSG",
                        "Age": "Allages",
                    },
                    "BAP1": {
                        "Actionability": "SA-CSG",
                        "Age": "Allages",
                    }
                },
                "Gene",
                "SAMPLE_ID"
            ),
            "# [SAMPLE_ID] Column format of gene_predispositions:\n" +
            "# [SAMPLE_ID]    a list of character strings of the form:\n" +
            "# [SAMPLE_ID]    [Gene]_[Actionability]_[Age]\n" +
            "# [SAMPLE_ID]    where:\n" +
            "# [SAMPLE_ID]        Actionability: [ HA-CSG | MA-CSG | SA-CSG ]\n" +
            "# [SAMPLE_ID]        Age: [ Age<30 | Allages ]\n"
        )
    ]
)
def test_generate_predisposition_column_format_description(inputs, want):
    # is there more testing needed here? like empty dict, or empty inner dict, or something?
    got = generate_predisposition_column_format_description(*inputs)
    assert got == want


def test_get_genomic_location():
    variant = "chromosome:position:ref>alt"
    out = "chromosome:position"
    assert get_genomic_location(variant) == out


def test_get_dna_change():
    variant = "chromosome:position:ref>alt"
    out = "ref>alt"
    assert get_dna_change(variant) == out


# def test_get_small_variant_calls_file_path()
# does this one need a test too?

# def test_lookup_predisposition_variants()


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

    source = "doi_url"

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
        f"# [{sample_id}]    file: {cancer_susceptibility_genes}\n"
        f"# [{sample_id}]    source: {source}\n"
        f"# [{sample_id}] Small variant calls are defined in: {small_variant_calls}\n"
        f"# [{sample_id}] Cumulative length of all the targeted coding regions (in millions of bases): {length_of_targeted_coding_regions:.2f}\n"
        f"# [{sample_id}] Tumor purity (as a fraction between 0 and 1): {tumor_purity:.2f}\n"
        f"# [{sample_id}]    The tumor purity is provided as an input parameter for tsoppy.\n"
        f"# [{sample_id}] Gene_predisposition column format:\n"
        f"# [{sample_id}]    [Gene]_[Actionability]_[Age]\n"
        f"# [{sample_id}]    where:\n"
        f"# [{sample_id}]        Actionability: [ HA-CSG | MA-CSG | SA-CSG ]\n"
        f"# [{sample_id}]        Age: [ Allages | Age<30 ]\n"
    )

    predisposition_column_format_description = (
        f"# [{sample_id}] Gene_predisposition column format:\n"
        f"# [{sample_id}]    [Gene]_[Actionability]_[Age]\n"
        f"# [{sample_id}]    where:\n"
        f"# [{sample_id}]        Actionability: [ HA-CSG | MA-CSG | SA-CSG ]\n"
        f"# [{sample_id}]        Age: [ Allages | Age<30 ]\n"
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
        source,
        cancer_susceptibility_genes,
        small_variant_calls,
        length_of_targeted_coding_regions,
        tumor_purity,
        pl.DataFrame(predisposition_variants),
        predisposition_column_format_description,
        output_file,
    )

    assert Path(expected_output_file).read_text() == Path(
        output_file).read_text()


# def test_generate_report()
