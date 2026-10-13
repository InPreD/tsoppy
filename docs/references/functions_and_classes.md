# Functions and classes

## Functions

### `tsoppy.general.file_parser.Parse_section_tsv`

The function parses a sectioned tsv file and returns potential headers and a dictionary containing the tables within the sections.

#### Input

1. `path: str`: Path to sectioned tsv file. The file may contain headers on top of the first section and the sections need to be separated by at least one empty line.

    ```tsv
    some header information
    more header information

    [section1]
    column1 column2
    1   one
    2   two
    3   three

    [section2]
    key1 value1
    key2 value2
    ```

1. `key_value_sections: list[str]`: List of strings representing sections that are known to be key value pairs. May be a list of length 0 if none of the section should be treated as key value pairs.

    ```python
    ['section2']
    ```

#### Output

1. `list[str]`: Each field on top of the first section is parsed if it is not an empty string, added to the list and returned. These fields usually contain header information.

    ```python
    ['some header information', 'more header information']
    ```
1. `dict[str, polars.DataFrame]`: A dictionary where the keys are the section names and the values are the tables (or any defined key value pairs) inside those sections represented as [data frames](https://docs.pola.rs/api/python/stable/reference/dataframe/index.html). The first row of each table is treated as column names unless the section is a list of key value pairs. Then the first column is converted to column names.

    ```python
    {
        'section1': polars.DataFrame({
            'column1': [1, 2, 3],
            'column2': ['one', 'two', 'three']
        }),
        'section2': polars.DataFrame({
            'key1': ['value1'],
            'key2': ['value2']
        })
    }
    ```
