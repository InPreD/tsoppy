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


test_data_dir = "tests/test_data/report_predispositions"


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
