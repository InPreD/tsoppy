# Development

If you are contributing to the codebase of tsoppy, we would like for you to read this section to familiarize yourself with the structure and guidelines of the project.

## Preparation

### Development environment setup

We recommend to install and use [**Visual Studio Code**](https://code.visualstudio.com/) for development but there might be other options more suitable if you are already familiar with them and if those support devcontainers. Download the appropriate version of Visual Studio Code and install it on your machine.

In order for you to use the devcontainer in this repository, you need to have [**Docker**](https://docs.docker.com/get-started/get-docker/) installed. Depending on your OS you should either go with Docker Desktop or install [Docker Engine](https://docs.docker.com/engine/install/) for linux distributions. Also, consider the [VS Code docs for setting up Docker](https://code.visualstudio.com/docs/devcontainers/containers#_installation).

### Clone repository

Open VS Code, open a new terminal, navigate to to directory you want to clone the repository into and run `git clone https://github.com/InPreD/tsoppy.git`. Now open the folder you just cloned to enter the local copy of the repository.

### Start devcontainer

Open the [`Command Palette`](https://code.visualstudio.com/api/ux-guidelines/command-palette) and choose `Dev Containers: Rebuild Container`. VS Code will launch the devcontainer which contains all the dependencies you need for developing in this repository. You can add extensions according to your needs and suggest changes to the devcontainer definition if they are relevant to other developers.

## Recommended workflow

1. Create an [issue](https://github.com/InPreD/tsoppy/issues/new) which is assigned to yourself, labeled as `enhancement` and of type `feature`. Include an adequate description and link to any relevant old TSOPPI code or other documents.
1. In the issue, use the `Create a branch` shortcut to create a new feature branch which should **always** branch off `develop`.
1. Check out your new branch locally and follow this guide for designing a subpackage.
1. Commit your changes often and as logically structured parts. Use [commit message conventions](https://inpred.github.io/24-03_bioinfo_ws/#19).
1. Before pushing, ensure that `ruff format --check` is happy with formatting.
1. When you are done with the work or want to get feedback, open a [pull request](https://github.com/InPreD/tsoppy/pulls). Your branches should be `base: develop` and `compare: <your feature branch>`. Assign yourself as assignee and include reviewers that can give you feedback on code quality and functionality.
1. Keep an eye on the github actions passing to make sure your unit tests work and your code is linted.
1. Address all comments from the reviewers and when your changes were approved it is time to merge.

## Recommended packages

The following packages are recommended to be used for the indicated purposes. We are open to discuss any of the following decisions depending on cost and benefit of the switch. If new purposes arise we will expand the list. Prior to suggesting new packages, please do some research on what people recommend and if the suggested package is actually maintained and fairly adopted by the community.

purpose | package
---|---
cli | [typer](https://typer.tiangolo.com/)
linting | [isort](https://pycqa.github.io/isort/), [ruff](https://docs.astral.sh/ruff/)
logging | [logging](https://docs.python.org/3/library/logging.html)
plotting | [plotnine](https://plotnine.org/)
read csv/tsv tables | [polars](https://docs.pola.rs/)
read json/toml/yaml | [msgspec](https://jcristharif.com/msgspec/index.html)
testing | [pytest](https://docs.pytest.org/en/stable/)
vcf | [CyVCF2](https://brentp.github.io/cyvcf2/)

## Repository structure

```txt
.
├── .devcontainer/ -> contains instructions to build and start a development container (https://containers.dev/)
│   ├── .devontainer.json
│   └── Dockerfile
├── .github -> github action scripts to perform CI related tasks such as testing, linting and building
├── .gitignore -> list of files to be excluded from git
├── docs/
│   ├── guides/
│   │   └── development.md -> this guide on how to contribute to tsoppy
│   ├── decisions/
│   │   └── <adr document name>.md -> architectual decision record
│   └── references/
│       ├── functions_and_classes.md -> description of how important functions and classes can be used
│       └── publications.md -> references to publications of used bioinfo packages, data sources, etc.
├── src/
│   └── tsoppy/ -> main package for tsoppy containing code for cli
│       ├── __init__.py
│       ├── cli.py
│       └── <subpackage>/ -> package containing a subfunctionality performing a specific task, typically linked to a subcommand
│           ├── __init__.py
│           └── <module>.py
├── tests/ -> unit tests and test data
│   ├── test_data/
│   │   └── <subpackage>_<module>/
│   │       └── <function>/
│   │           └── <test case name>.py
│   └── test_<subpackage>_<module>.py
├── CODEOWNERS -> specifying who is responsible
├── Dockerfile -> build recipe for docker
├── LICENSE -> license that we agreed on in InPreD bioinfo group
├── pyproject.toml -> project configuration, dependencies etc.
└── README.md -> documentation

```

## Designing a subcommand

1. Start by adding a folder to `src/tsoppy` using [snake_case](https://stringcase.org/cases/snake/) to name it. The name should be descriptive but also as short as possible. Any general functions and classes should be added to the existing subpackage `general`.
1. Initialize the folder with a `__init__.py` file.
1. If your subpackage consists of different parts create a module, a `.py` file named using snake_case, for each part. Otherwise, if only a single module is needed just call it `main.py`.
1. Include a short description on top of the file.
1. Any imports of python packages should follow and be sorted with `isort`. Also, add packages that are not installed yet to the `pyproject.toml`, section `[project]`, key `dependencies`.
1. Create a logger below the import section like so: `logger = logging.getLogger(__name__)`.
1. Functions should be named with snake_case. Starting with a capital letter indicates that the function is designed to be used outside of the module, e.g. `Public_function_name(...)`, while prefixing with `_` is for internal helper functions, e.g. `_internal_function(...)`. All input variables and output should be [typed](https://docs.python.org/3/library/typing.html). Include a short description under the function definition line. Also provide comments to describe individual sections of the function and for parts that are generally more complex.
1. Classes should be named with [PascalCase](https://stringcase.org/cases/pascal/). The class should contain a description and a list over all attributes. It also needs a constructor method `__init__` and potentially a [`__eq__` method](#alternatives-for-assert-value1--value2) for testing. The same guidelines listed for functions should be applied to methods.
1. Import your subpackage to `src/tsoppy/cli.py` like so: `from tsoppy.<subpackage>.<module> import <class or function>`. Subsequently, connect your subpackage to a command like so:

    ```python
    @app.command()
    def <subpackage name>(...):
        # call your function or create your class instance here
    ```

1. Remember to add [unit tests](#unit-testing).
1. In general, make things configurable and avoid hard-coding paths and variables that might be subjected to changes.

## Unit testing

We use pytest for unit testing and some kind of table-driven testing in order to reduce biolerplate code. Unit tests should be placed in `tests/<subpackage>_<module>_test.py`. If test data is necessary add it under `tests/test_data/<subpackage>_<module>/<function>/<test case name>.py`. Cover all edge cases as well as use cases from the different nodes.

Please find an example of a unit test below:

```python
from contextlib import nullcontext # in case of exception testing
from os import path # in case test data is used

from pytest import mark, raises # mark is required to use parametrize making the tests table-driven; raises is for exception testing

from tsoppy.<subpackage>.<module> import (
    <function>
    <class>
) # import the function(s) or class(es) you want to test

# Define path to test data - cannot be absolute due to different paths locally and in CI
test_data_dir = "tests/test_data/<subpackage>_<module>" # only required if test data is used; just to avoid repeating the path

@mark.parametrize(
    "inputs, exception, want", # exception can be skipped if no exceptions are tested
    [
        (
            # here should be a short description of the test case
            (
                "input1",
                path.join(test_data_dir, "<function>/<test_case>.<txt,tsv,csv,json,vcf>")
            ),
            nullcontext(),
            (
                "output1",
                2.0
            )
        ),
        ...
    ]
)
def test_<function>(inputs, exception, want):
    with exception: # if an exception is expected for some cases
        got = <function>(inputs[0], inputs[1])
        # assert that two values are equal is the simplest way to check if the output is as expected but not always adequate
        assert got[0] == want[0]
        assert got[1] == want[1]
```

> [!WARNING]
> Please keep in mind that the example is simplified and you need to adapt the tests to your needs. Start by defining your test cases and expected output and then modify the unit test accordingly.

### Potential exceptions to test for

- `FileNotFoundError`: If your function takes a file path as input.
- `polars.exceptions.NoDataError`: If you are reading a table into a polars dataframe.

### Alternatives for `assert value1 == value2`

- `dataframe1.equals(dataframe2)`: If you are comparing two polars dataframes.
- `__eq__(self, other)`: For classes, we can define a method to enable the `==` comparison:

    ```python
    def __eq__(self, other):
        if not isinstance(other, <class>):
            return False
        if self.<attribute1> != other.<attribute1>:
            return False
        ...
        return self.<attributeN> == other.<attributeN>:
    ```

## Running `tsoppy` cli

`tsoppy` needs to be reinstalled everytime you make changes to the code base which you want to test:

```bash
$ pip install .
```

Here are useful commands if you want to try out `tsoppy`:

```bash
$ tsoppy --help # help text for main command
$ tsoppy version # return version string
$ tsoppy <subcommand> --help # help text for sub command
```