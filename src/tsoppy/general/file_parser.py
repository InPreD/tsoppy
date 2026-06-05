import logging
import re

import polars

# Use logger that was set up in CLI
logger = logging.getLogger(__name__)


def Parse_section_tsv(
    path: str, key_value_sections: list[str]
) -> tuple[list[str], dict[str, polars.DataFrame]]:
    """Parse a sectioned TSV file into headers and a mapping of section names to DataFrames."""
    try:
        df = polars.read_csv(path, separator="\t", has_header=False)
    except FileNotFoundError:
        logger.error(f"File {path} not found.")
        raise
    except polars.exceptions.NoDataError:
        logger.error(f"File {path} is empty.")
        raise
    section_idx = _get_section_idx(df)
    headers = []
    if section_idx[0][1] != 1:
        headers = _parse_headers(df, section_idx[0][1] - 1)
    section_dfs = {}
    for section in section_idx:
        # create slice of dataframe for the section
        df_slice = df.slice(section[1], section[2])

        # check if row contains null values
        if any(item is None for item in df_slice.row(0)):
            df_slice = _handle_row_with_nulls(df_slice)

        # check if section is a key value section
        if section[0] in key_value_sections:
            # check that the section only contains two columns and transpose else log a warning
            if df_slice.width == 2:
                df_slice = df_slice.transpose()
            else:
                logger.warning(
                    f"Section {section[0]} is supposed to be a key value section but contains more than two columns."
                )

        # assume first row contains column names
        df_header = df_slice.head(1).to_dicts().pop()

        # remove first row and rename columns and link it to the section name
        section_dfs[section[0]] = df_slice.rename(df_header).slice(1)
    return headers, section_dfs


def _get_section_idx(df: polars.DataFrame) -> list[tuple[int, int]]:
    """Get the the name, start index and length of each section in the DataFrame."""
    section = ""
    section_start = 0
    section_length = 0
    section_idx = []
    for row in df.with_row_index().iter_rows():
        if row[1]:
            match = re.search(r"^\[(?P<section>.*)\]$", row[1])
            if match:
                section_start = row[0] + 1
                section = match.group("section")
        if all(item is None for item in row[1:]):
            section_length = row[0] - section_start
            if section_start > 0 and section_length > 0:
                section_idx.append((section, section_start, section_length))
                section_start = 0
                section_length = 0
        if row[0] == len(df) - 1:
            section_length = row[0] - section_start + 1
            if section_start > 0 and section_length > 0:
                section_idx.append((section, section_start, section_length))
    return section_idx


def _parse_headers(df: polars.DataFrame, header_rows: int) -> list[str]:
    """Parse the first rows from the top of the DataFrame as headers."""
    headers = []
    for row in df.head(header_rows).rows():
        for el in row:
            if el is not None:
                headers.append(el)
    return headers


def _handle_row_with_nulls(df: polars.DataFrame) -> polars.DataFrame:
    """Handle rows with null values removing empty columns and by filling with "-"."""

    # remove any columns that are completely null (no column header nor values)
    df = df.select(
        [polars.col(col)
         for col in df.columns if not df[col].null_count() == df.height]
    )

    # avoid null values by filling with "-"
    df = df.with_columns(polars.all().cast(polars.String).fill_null("-"))
    return df
