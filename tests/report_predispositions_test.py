"""
Report predispositions subpackage unit test.
"""

from contextlib import nullcontext
from pathlib import Path
from polars.exceptions import ColumnNotFoundError
from polars.testing import assert_frame_equal
from os import path

import polars as pl
import pytest

from tsoppy.report_predispositions.main import (
    check_columns_presence_in_file,
    check_gene_name_column_presence_in_column_list,
    format_predisposition_variants_before_printing,
    generate_predisposition_column_format_description,
    generate_report,
    get_small_variant_calls_file_path,
    get_workflow_id_string,
    load_data_from_cancer_susceptibility_genes_file,
    lookup_predisposition_variants,
    parse_metadata_to_get_source,
    print_metadata_to_output_file,
    print_output_file,
    validate_tumor_purity_range,
    validate_uniqueness,
    jinja
)

# Define path to test data - cannot be absolute due to different paths locally and in CI
test_data_dir = "tests/test_data/report_predispositions"


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
    # inputs: column_list, filename
    [
        (
            (
                [
                    "Gene",
                    "Actionability",
                    "Age",
                    "X"
                ],
                path.join(test_data_dir,
                          "check_columns_presence_in_file_missing_column.csv"),
            ),
            pytest.raises(ColumnNotFoundError),
            None
        ),
        (
            (
                [
                    "Gene",
                    "Actionability",
                    "Age"
                ],
                path.join(test_data_dir,
                          "check_columns_presence_in_file_expected_output.csv"),
            ),
            nullcontext(),
            None
        ),
    ]
)
def test_check_columns_presence_in_file(inputs, exception, want):
    with exception:
        got = check_columns_presence_in_file(*inputs)
        assert got == want


@pytest.mark.parametrize(
    "inputs, want",
    # inputs: filename, column_list, gene_name_column
    [
        (
            (
                Path(path.join(
                    test_data_dir, "load_data_from_cancer_susceptibility_genes_file_input.csv")),
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
def test_load_data_from_cancer_susceptibility_genes_file(inputs, want):
    got = load_data_from_cancer_susceptibility_genes_file(*inputs)
    assert got == want


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


@pytest.mark.parametrize(
    "inputs, exception, want",
    [
        (
            (
                "sample1",
                "config.yaml",
                path.join(test_data_dir, "nomenclature.yaml"),
                path.join(test_data_dir, "dragen/standard"),
            ),
            nullcontext(),
            Path(path.join(test_data_dir, "dragen/standard", "Logs_Intermediates",
                           "DnaDragenCaller", "sample1", "sample1.hard-filtered.gvcf.gz")),
        ),
        (
            (
                "sample1",
                "config.yaml",
                path.join(test_data_dir, "nomenclature.yaml"),
                path.join(test_data_dir, "localapp/standard"),
            ),
            nullcontext(),
            Path(path.join(test_data_dir, "localapp/standard", "Logs_Intermediates",
                           "VariantMatching", "sample1", "sample1_MergedSmallVariants.genome.vcf")),
        ),
        (
            (
                "sample1",
                "config.yaml",
                path.join(test_data_dir, "nomenclature.yaml"),
                path.join(test_data_dir,
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
                path.join(test_data_dir, "nomenclature.yaml"),
                path.join(test_data_dir, "dragen/standard"),
            ),
            nullcontext(),
            "dragen_2.6.2.4",
        ),
        (
            (
                "config.yaml",
                path.join(test_data_dir, "nomenclature.yaml"),
                path.join(test_data_dir, "localapp/standard"),
            ),
            nullcontext(),
            "localapp_ruo-2.2.0.12",
        ),
        (
            (
                "config.yaml",
                path.join(test_data_dir, "nomenclature.yaml"),
                path.join(test_data_dir,
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


@pytest.mark.parametrize(
    "inputs, want",
    # inputs:
    #   cancer_susceptibility_genes_dict,
    #   sample_id,
    #   config_yaml,
    #   inpred_nomenclature,
    #   root_path,
    # want: pl.DataFrame
    [
        (
            # no small variants overlapping the cances susceptibility genes
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
                "sample2",
                "config.yaml",
                path.join(test_data_dir, "nomenclature.yaml"),
                path.join(test_data_dir, "dragen/standard"),
            ),
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
            )
        ),
        (
            # one cancer susceptibility gene overlapping a small variant
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
                "sample3",
                "config.yaml",
                path.join(test_data_dir, "nomenclature.yaml"),
                path.join(test_data_dir, "dragen/standard"),
            ),
            pl.DataFrame(
                {
                    "variant_ID": ["chr1:2000000:A>G"],
                    "variant_type": ["SNV"],
                    "Illumina_variant_class": ["Germline_DB"],
                    "overlapping_genes":
                    [
                        ["TP53"]
                    ],
                    "AD":
                    [
                        [419, 19, 0]
                    ],
                    "DP": [438],
                    "AF":
                    [
                        [0.0434, 0]
                    ],
                    "gene_predispositions":
                    [
                        ["TP53_HA-CSG_Age<30"]
                    ],
                }
            )
        ),
        (
            # two cancer susceptibility genes overlapping a small variant
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
                "sample4",
                "config.yaml",
                path.join(test_data_dir, "nomenclature.yaml"),
                path.join(test_data_dir, "dragen/standard"),
            ),
            pl.DataFrame(
                {
                    "variant_ID": ["chr1:2000000:A>G"],
                    "variant_type": ["SNV"],
                    "Illumina_variant_class": ["Germline_DB"],
                    "overlapping_genes":
                    [
                        ["BRCA1", "TP53"]
                    ],
                    "AD":
                    [
                        [419, 19, 0]
                    ],
                    "DP": [438],
                    "AF":
                    [
                        [0.0434, 0]
                    ],
                    "gene_predispositions":
                    [
                        ["BRCA1_MA-CSG_Allages", "TP53_HA-CSG_Age<30"]
                    ],
                }
            )
        ),
        (
            # two cancer susceptibility genes overlapping a small variant and one other gene overlapping the variant too
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
                "sample5",
                "config.yaml",
                path.join(test_data_dir, "nomenclature.yaml"),
                path.join(test_data_dir, "dragen/standard"),
            ),
            pl.DataFrame(
                {
                    "variant_ID": ["chr1:2000000:A>G"],
                    "variant_type": ["SNV"],
                    "Illumina_variant_class": ["Germline_DB"],
                    "overlapping_genes":
                    [
                        ["BRCA1", "TP53", "XXX"]
                    ],
                    "AD":
                    [
                        [419, 19, 0]
                    ],
                    "DP": [438],
                    "AF":
                    [
                        [0.0434, 0]
                    ],
                    "gene_predispositions":
                    [
                        ["BRCA1_MA-CSG_Allages", "TP53_HA-CSG_Age<30"]
                    ],
                }
            )
        ),
    ]
)
def test_lookup_predisposition_variants(inputs, want):
    got = lookup_predisposition_variants(*inputs)
    assert_frame_equal(left=got, right=want, check_dtypes=False)


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
    #   predisposition_column_format_description,
    #   output_file
    # want:   filename
    [
        (
            (
                "sample_id",
                "tsoppy version 0.0.1",
                "dragen_2.6.2.4",
                "doi_url",
                Path("cancer_susceptibility_genes.tsv"),
                Path("small_variant_calls.vcf"),
                1.283,
                0.7,
                "# [sample_id] Column format of gene_predispositions:\n" +
                "# [sample_id]    a list of character strings of the form:\n" +
                "# [sample_id]    [Gene]_[Actionability]_[Age]\n" +
                "# [sample_id]    where:\n" +
                "# [sample_id]        Actionability: [ HA-CSG | MA-CSG | SA-CSG ]\n" +
                "# [sample_id]        Age: [ Age<30 | Allages ]\n",
                Path("output.tsv")
            ),
            Path(path.join(test_data_dir,
                 "print_metadata_to_output_file_expected_output.tsv")),
        ),


    ],
)
def test_print_metadata_to_output_file(tmp_path, inputs, want):

    # create temporary directory
    directory = tmp_path / "test_reportPredispositions_main"
    directory.mkdir(parents=True, exist_ok=True)
    out_path = directory / inputs[9]

    print_metadata_to_output_file(*inputs[:9], Path(out_path))
    assert Path(out_path).read_text() == Path(want).read_text()

    # remove the temporary file
    out_path.unlink()
    # remove the temp directory
    directory.rmdir()


@pytest.mark.parametrize(
    "inputs, want",
    # inputs:
    #   sample_id,
    #   pl.DataFrame(predisposition_variants),
    # want:
    #   pl.DataFrame(formated_predisposition_variants)
    [
        (
            (
                "sample1",
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
            ),
            pl.DataFrame(
                {
                    "sample_id": [],
                    "overlapping_genes": [],
                    "genomic_location": [],
                    "DNA_change": [],
                    "depth_tumor_DNA": [],
                    "AF_tumor_DNA": [],
                    "gene_predispositions": [],
                }
            ),
        ),
        (
            (
                "sample2",
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
            ),
            pl.DataFrame(
                {
                    "sample_id":
                        [
                            "sample2"
                        ],
                    "overlapping_genes":
                        [
                            "TP53"
                        ],
                    "genomic_location":
                        [
                            "chr1:1000000"
                        ],
                    "DNA_change":
                        [
                            "A>G"
                        ],
                    "depth_tumor_DNA":
                        [
                            190
                        ],
                    "AF_tumor_DNA":
                        [
                            0.0054
                        ],
                    "gene_predispositions":
                        [
                            "TP53_HA-CSG_Age<30",
                        ],
                }
            ),
        ),
        (
            (
                "sample3",
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
            ),
            pl.DataFrame(
                {
                    "sample_id":
                        [
                            "sample3"
                        ],
                    "overlapping_genes":
                        [
                            "TP53, BRCA1"
                        ],
                    "genomic_location":
                        [
                            "chr1:1000000"
                        ],
                    "DNA_change":
                        [
                            "A>G"
                        ],
                    "depth_tumor_DNA":
                        [
                            190
                        ],
                    "AF_tumor_DNA":
                        [
                            0.0054
                        ],
                    "gene_predispositions":
                        [
                            "TP53_HA-CSG_Age<30, BRCA1_MA-CSG_Allages"
                        ],
                }
            ),
        ),
        (
            (
                "sample4",
                # a variant overlapping with two predisposition genes and one gene which is not in the cancer susceptibility genes
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
                                ["TP53", "BRCA1", "XXX"],
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
            ),
            pl.DataFrame(
                {
                    "sample_id":
                        [
                            "sample4"
                        ],
                    "overlapping_genes":
                        [
                            "TP53, BRCA1, XXX"
                        ],
                    "genomic_location":
                        [
                            "chr1:1000000"
                        ],
                    "DNA_change":
                        [
                            "A>G"
                        ],
                    "depth_tumor_DNA":
                        [
                            190
                        ],
                    "AF_tumor_DNA":
                        [
                            0.0054
                        ],
                    "gene_predispositions":
                        [
                            "TP53_HA-CSG_Age<30, BRCA1_MA-CSG_Allages"
                        ],
                }
            ),
        ),
    ],
)
def test_format_predisposition_variants_before_printing(inputs, want):
    got = format_predisposition_variants_before_printing(*inputs)
    assert_frame_equal(left=got, right=want, check_dtypes=False)


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
    # want:   filepath
    [
        (
            (
                "SAMPLE_ID",
                "tsoppy version 0.0.1",
                "localapp_ruo-2.2.0.12",
                "doi_url",
                Path("cancer_susceptibility_genes.tsv"),
                Path("small_variant_calls.vcf"),
                1.283,
                0.7,
                # a variant overlapping with two predisposition genes and one gene which is not in the cancer susceptibility genes
                # the same as merged_table dataframe plus a column for the predisposition gene info
                # localapp pl.Dataframe input
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
                                ["TP53", "BRCA1", "XXX"],
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
                                [0.0054, 0.0],
                            ],
                        "gene_predispositions":
                            [
                                ["TP53_HA-CSG_Age<30", "BRCA1_MA-CSG_Allages"],
                            ],
                    }
                ),
                "# [SAMPLE_ID] Column format of gene_predispositions:\n" +
                "# [SAMPLE_ID]    a list of character strings of the form:\n" +
                "# [SAMPLE_ID]    [Gene]_[Actionability]_[Age]\n" +
                "# [SAMPLE_ID]    where:\n" +
                "# [SAMPLE_ID]        Actionability: [ HA-CSG | MA-CSG | SA-CSG ]\n" +
                "# [SAMPLE_ID]        Age: [ Age<30 | Allages ]\n",
                Path("output.tsv")
            ),
            Path(path.join(test_data_dir, "print_output_file_expected_output.tsv")),
        ),
    ],
)
def test_print_output_file(tmp_path, inputs, want):
    # create temporary directory
    directory = tmp_path / "test_reportPredispositions_main" / "test_print_output"
    directory.mkdir(parents=True, exist_ok=True)
    out_path = directory / inputs[10]
    print_output_file(*inputs[:10], out_path)
    assert out_path.read_text() == want.read_text()
    # remove the temporary file
    out_path.unlink()
    # remove the temp directory
    directory.rmdir()


@pytest.mark.parametrize(
    "inputs, want",
    # inputs:
    #   output_file,
    #   tsoppy_version,
    #   sample_id,
    #   config_yaml,
    #   inpred_nomenclature,
    #   root_path,
    #   tumor_purity,
    #   cancer_susceptibility_genes,
    #   csg_column_list,
    #   gene_name_column
    # want: filename
    [
        (
            (
                "output.tsv",
                "tsoppy version 0.0.1",
                "sample5",
                "config.yaml",
                path.join(test_data_dir, "nomenclature.yaml"),
                path.join(test_data_dir, "dragen/standard"),
                0.7,
                path.join(
                    test_data_dir, "cancer_susceptibility_genes.csv"),
                ["Gene", "Actionability", "Age"],
                "Gene"
            ),
            Path(path.join(test_data_dir, "generate_report_expected_output.tsv"))
        )
    ]
)
def test_generate_report(tmp_path, inputs, want):
    # create temporary directory
    directory = tmp_path / "test_reportPredispositions_main" / "test_generate_report"
    directory.mkdir(parents=True, exist_ok=True)
    out_path = directory / inputs[0]
    generate_report(out_path, *inputs[1:])
    assert out_path.read_text() == want.read_text()
    # remove the temporary file
    out_path.unlink()
    # remove the temp directory
    directory.rmdir()


@pytest.mark.parametrize(
    "inputs, want",
    [
        (
            (
                "sample_id",
                "tsoppy version 0.0.1",
                "dragen_2.6.2.4",
                Path("cancer_susceptibility_genes.tsv"),
                "doi_url",
                Path("small_variant_calls.vcf"),
                1.283,
                0.70,
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
                pl.DataFrame(
                    {
                        "sample_id":
                        [
                            "sample_id"
                        ],
                        "overlapping_genes":
                        [
                            "TP53, BRCA1, XXX"
                        ],
                        "genomic_location":
                        [
                            "chr1:1000000"
                        ],
                        "DNA_change":
                        [
                            "A>G"
                        ],
                        "depth_tumor_DNA":
                        [
                            190
                        ],
                        "AF_tumor_DNA":
                        [
                            0.0054
                        ],
                        "gene_predispositions":
                        [
                            "TP53_HA-CSG_Age<30, BRCA1_MA-CSG_Allages"
                        ],
                    }
                ),
            ),
            Path(path.join(test_data_dir,
                 "jinja_expected_output.tsv")),
        ),
    ]
)
def test_jinja(inputs, want):
    assert jinja(*inputs) == want.read_text()
