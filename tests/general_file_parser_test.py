import polars
from os import path
from pytest import mark, raises
from contextlib import nullcontext


from tsoppy.general.file_parser import Parse_section_tsv

# Define path to test data - cannot be absolute due to different paths locally and in CI
test_data_dir = "tests/test_data/general"


@mark.parametrize(
    "inputs, exception, want",
    [
        (
            (
                path.join(test_data_dir, "parse_section_tsv/standard.tsv"),
                []
            ),
            nullcontext(),
            (
                [
                    "header1"
                ],
                {
                    "section1": polars.DataFrame(
                        {
                            "col1": ["value1", "value2"],
                            "col2": ["value3", "value4"]
                        }
                    )
                }
            )
        ),
        (
            (
                path.join(test_data_dir,
                          "parse_section_tsv/multiple_sections.tsv"),
                []
            ),
            nullcontext(),
            (
                [
                    "header1"
                ],
                {
                    "section1": polars.DataFrame(
                        {
                            "col1": ["value1", "value2"],
                            "col2": ["value3", "value4"]
                        }
                    ),
                    "section2": polars.DataFrame(
                        {
                            "col1": ["value1", "value2"],
                            "col2": ["value3", "value4"]
                        }
                    )
                }
            )
        ),
        (
            (
                path.join(test_data_dir,
                          "parse_section_tsv/no_headers.tsv"),
                []
            ),
            nullcontext(),
            (
                [],
                {
                    "section1": polars.DataFrame(
                        {
                            "col1": ["value1", "value2"],
                            "col2": ["value3", "value4"]
                        }
                    )
                }
            )
        ),
        (
            (
                path.join(test_data_dir,
                          "parse_section_tsv/extra_empty_lines.tsv"),
                []
            ),
            nullcontext(),
            (
                [
                    "header1"
                ],
                {
                    "section1": polars.DataFrame(
                        {
                            "col1": ["value1", "value2"],
                            "col2": ["value3", "value4"]
                        }
                    )
                }
            )
        ),
        (
            (
                path.join(test_data_dir,
                          "parse_section_tsv/null_columns.tsv"),
                []
            ),
            nullcontext(),
            (
                [
                    "header1"
                ],
                {
                    "section1": polars.DataFrame(
                        {
                            "col1": ["value1", "value2"],
                            "col2": ["value3", "value4"]
                        }
                    )
                }
            )
        ),
        (
            (
                path.join(test_data_dir,
                          "parse_section_tsv/empty_first_column_name.tsv"),
                []
            ),
            nullcontext(),
            (
                [
                    "header1"
                ],
                {
                    "section1": polars.DataFrame(
                        {
                            "-": ["value1", "value2"],
                            "col2": ["value3", "value4"]
                        }
                    )
                }
            )
        ),
        (
            (
                path.join(test_data_dir,
                          "parse_section_tsv/key_value.tsv"),
                [
                    "section1"
                ]
            ),
            nullcontext(),
            (
                [
                    "header1"
                ],
                {
                    "section1": polars.DataFrame(
                        {
                            "key1": ["value1"],
                            "key2": ["value2"]
                        }
                    )
                }
            )
        ),
        (
            (
                path.join(test_data_dir,
                          "parse_section_tsv/non-existent.tsv"),
                []
            ),
            raises(FileNotFoundError),
            (
                [],
                {}
            )
        ),
        (
            (
                path.join(test_data_dir,
                          "parse_section_tsv/empty.tsv"),
                []
            ),
            raises(polars.exceptions.NoDataError),
            (
                [],
                {}
            )
        )
    ]
)
def test_parse_section_tsv(inputs, exception, want):
    with exception:
        got = Parse_section_tsv(inputs[0], inputs[1])
        assert got[0] == want[0]
        for key in want[1].keys():
            assert key in got[1]
            assert got[1][key].equals(want[1][key])
