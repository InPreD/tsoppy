"""
Report predispositions subpackage unit test.
"""

from contextlib import nullcontext
from pathlib import Path
from polars.exceptions import (ColumnNotFoundError, NoDataError)
from os import path

import importlib.metadata
import polars as pl
import pytest

from tsoppy.general.classes import (
    SmallVariantGenomeVcf,
    WorkflowOutput,
)

from tsoppy.report_predispositions.main import (
    check_columns_presence_in_file,
    check_file_existence_and_content_presence,
    check_gene_name_column_presence_in_column_list,
    generate_predisposition_column_format_description,
    generate_report,
    get_small_variant_calls_file_path,
    get_workflow_id_string,
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


# Define path to test data - cannot be absolute due to different paths locally and in CI
test_data_dir_general_classes = "tests/test_data/general_classes"


@pytest.mark.parametrize(
    "inputs, exception, want",
    [
        (
            (
                "sample1",
                "config.yaml",
                path.join(test_data_dir_general_classes, "nomenclature.yaml"),
                path.join(test_data_dir_general_classes, "dragen/standard"),
            ),
            nullcontext(),
            Path(path.join(test_data_dir_general_classes, "dragen/standard", "Logs_Intermediates",
                           "DnaDragenCaller", "sample1", "sample1.hard-filtered.gvcf.gz")),
        ),
        (
            (
                "sample1",
                "config.yaml",
                path.join(test_data_dir_general_classes, "nomenclature.yaml"),
                path.join(test_data_dir_general_classes, "localapp/standard"),
            ),
            nullcontext(),
            Path(path.join(test_data_dir_general_classes, "localapp/standard", "Logs_Intermediates",
                           "VariantMatching", "sample1", "sample1_MergedSmallVariants.genome.vcf")),
        ),
        (
            (
                "sample1",
                "config.yaml",
                path.join(test_data_dir_general_classes, "nomenclature.yaml"),
                path.join(test_data_dir_general_classes,
                          "dragen/non-existent"),
            ),
            pytest.raises(FileNotFoundError),
            None,
        ),
    ],
)
def test_get_small_variant_calls_file_path(inputs, exception, want):
    with exception:
        got = get_small_variant_calls_file_path(*inputs)
        assert got == want


@pytest.mark.parametrize(
    "inputs, exception, want",
    [
        (
            (
                "config.yaml",
                path.join(test_data_dir_general_classes, "nomenclature.yaml"),
                path.join(test_data_dir_general_classes, "dragen/standard"),
            ),
            nullcontext(),
            "dragen_2.6.2.4",
        ),
        (
            (
                "config.yaml",
                path.join(test_data_dir_general_classes, "nomenclature.yaml"),
                path.join(test_data_dir_general_classes, "localapp/standard"),
            ),
            nullcontext(),
            "localapp_ruo-2.2.0.12",
        ),
        (
            (
                "config.yaml",
                path.join(test_data_dir_general_classes, "nomenclature.yaml"),
                path.join(test_data_dir_general_classes,
                          "dragen/non-existent"),
            ),
            pytest.raises(FileNotFoundError),
            None,
        ),
    ],
)
def test_get_workflow_id_string(inputs, exception, want):
    with exception:
        got = get_workflow_id_string(*inputs)
        assert got == want


# def test_lookup_predisposition_variants()


@pytest.mark.parametrize(
    "inputs, want",
    # inputs:
    #   sample_id,
    #   tsoppy_version,
    #   workflow_id,
    #   source,
    #   cancer_susceptibility_genes,
    #   small_variant_calls,
    #   coding_region_size_in_megabases,
    #   tumor_purity,
    #   pl.DataFrame(predisposition_variants),
    #   predisposition_column_format_description,
    #   output_file
    # want:   filename, metadata, header + variants
    [
        (
            (
                "sample1",
                "tsoppy version 0.0.1",
                "dragen_2.6.2.4",
                "doi_url",
                Path("cancer_susceptibility_genes.tsv"),
                Path("small_variant_calls.vcf"),
                1.283,
                0.7,
                # empty dataframe
                # the same as merged_table dataframe plus a column for the predisposition gene info
                pl.DataFrame(
                    {
                        "variant_ID": [],
                        "variant_type": [],
                        "Illumina_variant_class": [],
                        "overlapping_genes": [],
                        "AD": [],
                        "DP": [],
                        "AF": [],
                        "gene_predispositions": [],
                    }
                ),
                "# [sample1] Column format of gene_predispositions:\n" +
                "# [sample1]    a list of character strings of the form:\n" +
                "# [sample1]    [Gene]_[Actionability]_[Age]\n" +
                "# [sample1]    where:\n" +
                "# [sample1]        Actionability: [ HA-CSG | MA-CSG | SA-CSG ]\n" +
                "# [sample1]        Age: [ Age<30 | Allages ]\n",
                Path("output_sample1.tsv")
            ),
            (
                Path("expected_output_sample1.tsv"),
                (
                    f"# [sample1] Pipeline: tsoppy version 0.0.1\n"
                    f"# [sample1] Workflow id: dragen_2.6.2.4\n"
                    f"# [sample1] Cancer susceptibility genes are defined in:\n"
                    f"# [sample1]    file: cancer_susceptibility_genes.tsv\n"
                    f"# [sample1]    source: doi_url\n"
                    f"# [sample1] Small variant calls are defined in: small_variant_calls.vcf\n"
                    f"# [sample1] Cumulative length of all the targeted coding regions (in millions of bases, seq. depth >= 50): 1.28\n"
                    f"# [sample1] Tumor purity (as a fraction between 0 and 1): 0.70\n"
                    f"# [sample1]    The tumor purity is provided as an input parameter to tsoppy.\n"
                    f"# [sample1] Column format of gene_predispositions:\n"
                    f"# [sample1]    a list of character strings of the form:\n"
                    f"# [sample1]    [Gene]_[Actionability]_[Age]\n"
                    f"# [sample1]    where:\n"
                    f"# [sample1]        Actionability: [ HA-CSG | MA-CSG | SA-CSG ]\n"
                    f"# [sample1]        Age: [ Age<30 | Allages ]\n"
                ),
                (
                    f"variant_ID\t"
                    f"variant_type\t"
                    f"Illumina_variant_class\t"
                    f"overlapping_genes\t"
                    f"AD\t"
                    f"DP\t"
                    f"AF\t"
                    f"gene_predispositions\n"
                ),
            ),
        ),
        (
            (
                "sample2",
                "tsoppy version 0.0.1",
                "dragen_2.6.2.4",
                "doi_url",
                Path("cancer_susceptibility_genes.tsv"),
                Path("small_variant_calls.vcf"),
                1.283,
                0.7,
                # a variant overlapping with one predisposition gene
                # the same as merged_table dataframe plus a column for the predisposition gene info
                pl.DataFrame(
                    {
                        "variant_ID":
                            [
                                "chr1:1000000:A>G",
                            ],
                        "variant_type":
                            [
                                "SNV",
                            ],
                        "Illumina_variant_class":
                            [
                                "VCF_filtered",
                            ],
                        "overlapping_genes":
                            [
                                ["TP53"],
                            ],
                        "AD":
                            [
                                [189, 1, 0],
                            ],
                        "DP":
                            [
                                190,
                            ],
                        "AF":
                            [
                                [0.0054, 0.0]
                            ],
                        "gene_predispositions":
                            [
                                "TP53_HA-CSG_Age<30",
                            ],
                    }
                ),
                "# [sample2] Column format of gene_predispositions:\n" +
                "# [sample2]    a list of character strings of the form:\n" +
                "# [sample2]    [Gene]_[Actionability]_[Age]\n" +
                "# [sample2]    where:\n" +
                "# [sample2]        Actionability: [ HA-CSG | MA-CSG | SA-CSG ]\n" +
                "# [sample2]        Age: [ Age<30 | Allages ]\n",
                Path("output_sample2.tsv")
            ),
            (
                Path("expected_output_sample2.tsv"),
                (
                    f"# [sample2] Pipeline: tsoppy version 0.0.1\n"
                    f"# [sample2] Workflow id: dragen_2.6.2.4\n"
                    f"# [sample2] Cancer susceptibility genes are defined in:\n"
                    f"# [sample2]    file: cancer_susceptibility_genes.tsv\n"
                    f"# [sample2]    source: doi_url\n"
                    f"# [sample2] Small variant calls are defined in: small_variant_calls.vcf\n"
                    f"# [sample2] Cumulative length of all the targeted coding regions (in millions of bases, seq. depth >= 50): 1.28\n"
                    f"# [sample2] Tumor purity (as a fraction between 0 and 1): 0.70\n"
                    f"# [sample2]    The tumor purity is provided as an input parameter to tsoppy.\n"
                    f"# [sample2] Column format of gene_predispositions:\n"
                    f"# [sample2]    a list of character strings of the form:\n"
                    f"# [sample2]    [Gene]_[Actionability]_[Age]\n"
                    f"# [sample2]    where:\n"
                    f"# [sample2]        Actionability: [ HA-CSG | MA-CSG | SA-CSG ]\n"
                    f"# [sample2]        Age: [ Age<30 | Allages ]\n"
                ),
                (
                    f"variant_ID\t"
                    f"variant_type\t"
                    f"Illumina_variant_class\t"
                    f"overlapping_genes\t"
                    f"AD\t"
                    f"DP\t"
                    f"AF\t"
                    f"gene_predispositions\n"
                    f"chr1:1000000:A>G\t"
                    f"SNV\t"
                    f"VCF_filtered\t"
                    f"TP53\t"
                    f"189, 1, 0\t"
                    f"190\t"
                    f"0.0054, 0.0\t"
                    f"TP53_HA-CSG_Age<30\n"
                ),
            ),
        ),
        (
            (
                "sample3",
                "tsoppy version 0.0.1",
                "dragen_2.6.2.4",
                "doi_url",
                Path("cancer_susceptibility_genes.tsv"),
                Path("small_variant_calls.vcf"),
                1.283,
                0.7,
                # a variant overlapping with two predisposition genes
                # the same as merged_table dataframe plus a column for the predisposition gene info
                pl.DataFrame(
                    {
                        "variant_ID":
                            [
                                "chr1:1000000:A>G",
                            ],
                        "variant_type":
                            [
                                "SNV",
                            ],
                        "Illumina_variant_class":
                            [
                                "VCF_filtered",
                            ],
                        "overlapping_genes":
                            [
                                ["TP53", "BRCA1"],
                            ],
                        "AD":
                            [
                                [189, 1, 0],
                            ],
                        "DP":
                            [
                                190,
                            ],
                        "AF":
                            [
                                [0.0054, 0.0]
                            ],
                        "gene_predispositions":
                            [
                                ["TP53_HA-CSG_Age<30", "BRCA1_MA-CSG_Allages"],
                            ],
                    }
                ),
                "# [sample3] Column format of gene_predispositions:\n" +
                "# [sample3]    a list of character strings of the form:\n" +
                "# [sample3]    [Gene]_[Actionability]_[Age]\n" +
                "# [sample3]    where:\n" +
                "# [sample3]        Actionability: [ HA-CSG | MA-CSG | SA-CSG ]\n" +
                "# [sample3]        Age: [ Age<30 | Allages ]\n",
                Path("output_sample3.tsv")
            ),
            (
                Path("expected_output_sample3.tsv"),
                (
                    f"# [sample3] Pipeline: tsoppy version 0.0.1\n"
                    f"# [sample3] Workflow id: dragen_2.6.2.4\n"
                    f"# [sample3] Cancer susceptibility genes are defined in:\n"
                    f"# [sample3]    file: cancer_susceptibility_genes.tsv\n"
                    f"# [sample3]    source: doi_url\n"
                    f"# [sample3] Small variant calls are defined in: small_variant_calls.vcf\n"
                    f"# [sample3] Cumulative length of all the targeted coding regions (in millions of bases, seq. depth >= 50): 1.28\n"
                    f"# [sample3] Tumor purity (as a fraction between 0 and 1): 0.70\n"
                    f"# [sample3]    The tumor purity is provided as an input parameter to tsoppy.\n"
                    f"# [sample3] Column format of gene_predispositions:\n"
                    f"# [sample3]    a list of character strings of the form:\n"
                    f"# [sample3]    [Gene]_[Actionability]_[Age]\n"
                    f"# [sample3]    where:\n"
                    f"# [sample3]        Actionability: [ HA-CSG | MA-CSG | SA-CSG ]\n"
                    f"# [sample3]        Age: [ Age<30 | Allages ]\n"
                ),
                (
                    f"variant_ID\t"
                    f"variant_type\t"
                    f"Illumina_variant_class\t"
                    f"overlapping_genes\t"
                    f"AD\t"
                    f"DP\t"
                    f"AF\t"
                    f"gene_predispositions\n"
                    f"chr1:1000000:A>G\t"
                    f"SNV\t"
                    f"VCF_filtered\t"
                    f"TP53, BRCA1\t"
                    f"189, 1, 0\t"
                    f"190\t"
                    f"0.0054, 0.0\t"
                    f"TP53_HA-CSG_Age<30, BRCA1_MA-CSG_Allages\n"
                ),
            ),
        ),
    ],
)
def test_print_predisposition_variants_to_output_file(tmp_path, inputs, want):

    # create temporary directory
    directory = tmp_path / "test_reportPredispositions_main"
    directory.mkdir(parents=True, exist_ok=True)

    out_path = directory / inputs[10]
    print_predisposition_variants_to_output_file(*inputs[:10], out_path)

    expected_path = directory / want[0]
    expected_content = want[1] + want[2]
    Path(expected_path).write_text(expected_content)

    assert Path(out_path).read_text() == Path(expected_path).read_text()

    # remove the temporary files
    out_path.unlink()
    expected_path.unlink()

    # remove the temp directory
    directory.rmdir()


# def test_generate_report()
