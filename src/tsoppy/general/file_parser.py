import logging
import re

import polars

# Use logger that was set up in CLI
logger = logging.getLogger(__name__)


class sectionIdx:
    empty = False
    length = 0
    start = 0

    def __init__(self, name: str, start: int, length: int = 0):
        self.name = name
        self.start = start
        self.length = length

    def __eq__(self, other):
        if not isinstance(other, sectionIdx):
            return False
        attr_to_compare = ["empty", "length", "name", "start"]
        return all(vars(self).get(k) == vars(other).get(k) for k in attr_to_compare)

    def finalize(self):
        self.set_length(self.length)
        return self

    def included(self, idx_list: list[sectionIdx]) -> bool:
        # Ignore empty sections
        if self.name == "":
            return True
        return any(item.name == self.name for item in idx_list)

    def set_length(self, length: int):
        if length == 0:
            self.empty = True
        self.length = length


def Parse_section_tsv(
    path: str, key_value_sections: list[str]
) -> tuple[list[str], dict[str, polars.DataFrame]]:
    """Parse a sectioned TSV file into headers and a mapping of section names to DataFrames."""
    return _parse_section_sep_val(path, key_value_sections, "\t")


def Parse_section_csv(
    path: str, key_value_sections: list[str]
) -> tuple[list[str], dict[str, polars.DataFrame]]:
    """Parse a sectioned CSV file into headers and a mapping of section names to DataFrames."""
    return _parse_section_sep_val(path, key_value_sections, ",")


def _parse_section_sep_val(
    path: str, key_value_sections: list[str], sep: str
) -> tuple[list[str], dict[str, polars.DataFrame]]:
    """Parse a sectioned file containing separated values (TSV, CSV, etc.) into headers and a mapping of section names to DataFrames."""
    try:
        df = polars.read_csv(path, separator=sep, has_header=False)
    except FileNotFoundError:
        logger.error(f"File {path} not found.")
        raise
    except polars.exceptions.NoDataError:
        logger.error(f"File {path} is empty.")
        raise
    section_idx = _get_section_idx(df)
    headers = []
    if section_idx[0].start != 1:
        headers = _parse_headers(df, section_idx[0].start - 1)
    section_dfs = {}
    for section in section_idx:
        # create slice of dataframe for the section
        if section.empty:
            section_dfs[section.name] = polars.DataFrame()
            continue
        df_slice = df.slice(section.start, section.length)

        # check if row contains null values
        if any(item is None for item in df_slice.row(0)):
            df_slice = _handle_row_with_nulls(df_slice)

        # check if section is a key value section
        if section.name in key_value_sections:
            # check that the section only contains two columns and transpose else log a warning
            if df_slice.width == 2:
                df_slice = df_slice.transpose()
            else:
                logger.warning(
                    f"Section {section.name} is supposed to be a key value section but contains more than two columns."
                )

        # assume first row contains column names
        df_header = df_slice.head(1).to_dicts().pop()

        # remove first row and rename columns and link it to the section name
        section_dfs[section.name] = df_slice.rename(df_header).slice(1)
    return headers, section_dfs


def _get_section_idx(df: polars.DataFrame) -> list[tuple[sectionIdx]]:
    """Get the the name, start index and length of each section in the DataFrame."""
    section = sectionIdx("", 0)
    section_idx = []
    for row in df.with_row_index().iter_rows():
        if row[1]:
            match = re.search(r"^\[(?P<section>.*)\]$", row[1])
            if match:
                if not section.included(section_idx):
                    section_idx.append(section.finalize())
                section = sectionIdx(match.group("section"), row[0] + 1)
                continue
            else:
                section.set_length(row[0] - section.start + 1)
        if all(item is None for item in row[1:]):
            continue
    if not section.included(section_idx):
        section_idx.append(section.finalize())
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
    """Handle rows with null values by removing empty columns and replacing remaining missing values with "-"."""

    # remove any columns that are completely null (no column header nor values)
    df = df.select(
        [polars.col(col) for col in df.columns if not df[col].null_count() == df.height]
    )

    # avoid null values by filling with "-"
    df = df.with_columns(polars.all().cast(polars.String).fill_null("-"))
    return df
