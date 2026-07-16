from contextlib import nullcontext
from os import path

import polars
from pytest import mark, raises

from tsoppy.general.file_parser import (
    sectionIdx,
    _parse_section_sep_val,
    _get_section_idx,
    _handle_row_with_nulls,
    _parse_headers,
)

# Define path to test data - cannot be absolute due to different paths locally and in CI
test_data_dir = "tests/test_data/general_file_parser"


@mark.parametrize(
    "inputs, exception, want",
    [
        (
            # Standard tsv case with one section and headers
            (
                path.join(test_data_dir, "parse_section_sep_val/standard_tsv.tsv"),
                [],
                "\t",
            ),
            nullcontext(),
            (
                ["header1"],
                {
                    "section1": polars.DataFrame(
                        {"col1": ["value1", "value2"], "col2": ["value3", "value4"]}
                    )
                },
            ),
        ),
        (
            # Standard csv case with one section and headers
            (
                path.join(test_data_dir, "parse_section_sep_val/standard_csv.tsv"),
                [],
                ",",
            ),
            nullcontext(),
            (
                ["header1"],
                {
                    "section1": polars.DataFrame(
                        {"col1": ["value1", "value2"], "col2": ["value3", "value4"]}
                    )
                },
            ),
        ),
        (
            # Standard case with multiple sections and headers
            (
                path.join(test_data_dir, "parse_section_sep_val/multiple_sections.tsv"),
                [],
                "\t",
            ),
            nullcontext(),
            (
                ["header1"],
                {
                    "section1": polars.DataFrame(
                        {"col1": ["value1", "value2"], "col2": ["value3", "value4"]}
                    ),
                    "section2": polars.DataFrame(
                        {"col1": ["value1", "value2"], "col2": ["value3", "value4"]}
                    ),
                },
            ),
        ),
        (
            # No headers
            (
                path.join(test_data_dir, "parse_section_sep_val/no_headers.tsv"),
                [],
                "\t",
            ),
            nullcontext(),
            (
                [],
                {
                    "section1": polars.DataFrame(
                        {"col1": ["value1", "value2"], "col2": ["value3", "value4"]}
                    )
                },
            ),
        ),
        (
            # Empty section
            (
                path.join(test_data_dir, "parse_section_sep_val/empty_section.tsv"),
                [],
                "\t",
            ),
            nullcontext(),
            (
                [],
                {
                    "empty": polars.DataFrame(),
                    "section1": polars.DataFrame(
                        {"col1": ["value1", "value2"], "col2": ["value3", "value4"]}
                    ),
                },
            ),
        ),
        (
            # Extra empty lines between sections and headers
            (
                path.join(test_data_dir, "parse_section_sep_val/extra_empty_lines.tsv"),
                [],
                "\t",
            ),
            nullcontext(),
            (
                ["header1"],
                {
                    "section1": polars.DataFrame(
                        {"col1": ["value1", "value2"], "col2": ["value3", "value4"]}
                    )
                },
            ),
        ),
        (
            # Columns containing only null values
            (
                path.join(test_data_dir, "parse_section_sep_val/null_columns.tsv"),
                [],
                "\t",
            ),
            nullcontext(),
            (
                ["header1"],
                {
                    "section1": polars.DataFrame(
                        {"col1": ["value1", "value2"], "col2": ["value3", "value4"]}
                    )
                },
            ),
        ),
        (
            # Empty first column name
            (
                path.join(
                    test_data_dir, "parse_section_sep_val/empty_first_column_name.tsv"
                ),
                [],
                "\t",
            ),
            nullcontext(),
            (
                ["header1"],
                {
                    "section1": polars.DataFrame(
                        {"-": ["value1", "value2"], "col2": ["value3", "value4"]}
                    )
                },
            ),
        ),
        (
            # Key-value pairs instead of tabular data
            (
                path.join(test_data_dir, "parse_section_sep_val/key_value.tsv"),
                ["section1"],
                "\t",
            ),
            nullcontext(),
            (
                ["header1"],
                {
                    "section1": polars.DataFrame(
                        {"key1": ["value1"], "key2": ["value2"]}
                    )
                },
            ),
        ),
        (
            # Non-existent file
            (
                path.join(test_data_dir, "parse_section_sep_val/non-existent.tsv"),
                [],
                "\t",
            ),
            raises(FileNotFoundError),
            ([], {}),
        ),
        (
            # Empty file
            (path.join(test_data_dir, "parse_section_sep_val/empty.tsv"), [], "\t"),
            raises(polars.exceptions.NoDataError),
            ([], {}),
        ),
    ],
)
def test_parse_section_tsv(inputs, exception, want):
    with exception:
        got = _parse_section_sep_val(inputs[0], inputs[1], inputs[2])
        assert got[0] == want[0]
        for key in want[1].keys():
            assert key in got[1]
            assert got[1][key].equals(want[1][key])


@mark.parametrize(
    "input, want",
    [
        (
            # Standard case with one section and no headers
            polars.DataFrame(
                {
                    "col1": ["[section1]", "col1", "value1"],
                    "col2": [None, "col2", "value2"],
                }
            ),
            [sectionIdx("section1", 1, 2)],
        ),
        (
            # Standard case with one section and headers
            polars.DataFrame(
                {
                    "col1": ["header1", None, "[section1]", "col1", "value1"],
                    "col2": [None, None, None, "col2", "value2"],
                }
            ),
            [sectionIdx("section1", 3, 2)],
        ),
        (
            # Extra empty lines prior the top section, no headers
            polars.DataFrame(
                {
                    "col1": [None, None, "[section1]", "col1", "value1"],
                    "col2": [None, None, None, "col2", "value2"],
                }
            ),
            [sectionIdx("section1", 3, 2)],
        ),
        (
            # Missing column name in one section
            polars.DataFrame(
                {
                    "col1": [None, None, "[section1]", None, "value1"],
                    "col2": [None, None, None, "col2", "value2"],
                }
            ),
            [sectionIdx("section1", 3, 2)],
        ),
        (
            # Multiple sections with extra empty lines
            polars.DataFrame(
                {
                    "col1": [
                        None,
                        None,
                        "[section1]",
                        "col1",
                        "value1",
                        None,
                        "[section2]",
                        "col1",
                        "value1",
                    ],
                    "col2": [
                        None,
                        None,
                        None,
                        "col2",
                        "value2",
                        None,
                        None,
                        "col2",
                        "value2",
                    ],
                }
            ),
            [sectionIdx("section1", 3, 2), sectionIdx("section2", 7, 2)],
        ),
        (
            # Column only contains null values
            polars.DataFrame(
                {
                    "col1": ["[section1]", "col1", "value1"],
                    "col2": [None, "col2", "value2"],
                    "col3": [None, None, None],
                }
            ),
            [sectionIdx("section1", 1, 2)],
        ),
    ],
)
def test_get_section_idx(input, want):
    got = _get_section_idx(input)
    assert got == want


@mark.parametrize(
    "inputs, want",
    [
        (
            # Standard case with headers in the first rows
            (
                polars.DataFrame(
                    {"col1": [None, "header2"], "col2": ["header1", None]}
                ),
                2,
            ),
            ["header1", "header2"],
        )
    ],
)
def test_parse_headers(inputs, want):
    got = _parse_headers(inputs[0], inputs[1])
    assert got == want


@mark.parametrize(
    "input, want",
    [
        (
            # Column containing only null values
            polars.DataFrame(
                {
                    "col1": ["col1", "value1"],
                    "col2": ["col2", "value2"],
                    "col3": [None, None],
                }
            ),
            polars.DataFrame({"col1": ["col1", "value1"], "col2": ["col2", "value2"]}),
        ),
        (
            # Column name is null and column containing only null values
            polars.DataFrame(
                {
                    "col1": [None, "value1"],
                    "col2": ["col2", "value2"],
                    "col3": [None, None],
                }
            ),
            polars.DataFrame({"col1": ["-", "value1"], "col2": ["col2", "value2"]}),
        ),
    ],
)
def test_handle_row_with_nulls(input, want):
    got = _handle_row_with_nulls(input)
    assert got.equals(want)
