# tsoppy

TSO500 v2 post processing cli

## Development

If you are contributing to the codebase of tsoppy, we would like for you to read this section to familiarize yourself with the structure and guidelines of the project.

### Recommended workflow

1. Create an [issue](https://github.com/InPreD/tsoppy/issues/new) which is assigned to yourself, labeled as `enhancement` and of type `feature`. Include an adequate description and link to any relevant old TSOPPI code or other documents.
1. In the issue, use the `Create a branch` shortcut to create a new feature branch which should **always** branch off `develop`.
1. Check out your new branch locally and follow this guide for designing a subpackage.
1. Commit your changes often and as logically structured parts. Use [commit message conventions](https://inpred.github.io/24-03_bioinfo_ws/#19).
1. Before pushing, ensure that `ruff format --check` is happy with formatting.
1. When you are done with the work or want to get feedback, open a [pull request](https://github.com/InPreD/tsoppy/pulls). Your branches should be `base: develop` and `compare: <your feature branch>`. Assign yourself as assignee and include reviewers that can give you feedback on code quality and functionality.
1. Keep an eye on the github actions passing to make sure your unit tests work and your code is linted.
1. Address all comments from the reviewers and when your changes were approved it is time to merge.

### Recommended packages

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

### Repository structure

```txt
.
├── .devcontainer/ -> contains instructions to build and start a development container (https://containers.dev/)
│   ├── .devontainer.json
│   └── Dockerfile
├── .github -> github action scripts to perform CI related tasks such as testing, linting and building
├── .gitignore -> list of files to be excluded from git
├── docs/
│   ├── decisions/
│   │    <adr document name>.md -> architectual decision record
│   └── references/
│       └── publications.md -> references to publications of used bioinfo packages, data sources, etc.
├── src/
│   └── tsoppy/ -> main package for tsoppy containing code for cli
│       ├── __init__.py
│       ├── cli.py
│       └── <subpackage>/ -> package containing a subfunctionality performing a specific task, typically linked to a subcommand
│           ├── __init__.py
│           └── main.py
├── tests/ -> unit tests and test data
│   ├── test_data/
│   │   └── <subpackage>_module/
│   │       └── <test case name>.py
│   └── test_<subpackage>_<module>.py
├── CODEOWNERS -> specifying who is responsible
├── Dockerfile -> build recipe for docker
├── LICENSE -> license that we agreed on in InPreD bioinfo group
├── pyproject.toml -> project configuration, dependencies etc.
└── README.md -> documentation

```

### Designing a subcommand

1. Start by adding a folder to `src/tsoppy` using [snake_case](https://stringcase.org/cases/snake/) to name it. The name should be descriptive but also as short as possible.
1. Initialize the folder with a `__init__.py` file.
1. If your subpackage consists of different parts create a module, a `.py` file named using snake_case, for each part. Otherwise, if only a single module is needed just call it `main.py`.
1. Include a short description on top of the file.
1. Any imports of packages should follow and be sorted with `isort`. Also, add packages that are not installed yet to the `pyproject.toml`.
1. Create a logger below the import section like so: `logger = logging.getLogger(__name__)`.
1. Functions should be named with snake_case. All input variables and output should be [typed](https://docs.python.org/3/library/typing.html). Include a short description under the function definition line. Also provide comments to describe individual sections of the function and for parts that are generally more complex.
1. Classes should be named with [PascalCase](https://stringcase.org/cases/pascal/). The class should contain a description and a list over all attributes. It also needs a constructor method `__init__`. The same guidelines listed for functions should be applied to methods.
1. Import your subpackage to `src/tsoppy/cli.py` like so: `from tsoppy.<subpackage>.<module> import <class or function>`. Subsequently, connect your subpackage to a command like so:

    ```python
    @app.command()
    def <subpackage name>():
    ```

1. Provide unit tests in `tests/test_<subpackage>_<module>.py`. If test data is necessary add it under `tests/test_data/<subpackage>_<module>/<test case name>`. Cover all edge cases as well as use cases from the different nodes.
1. In general, make things configurable and avoid hard-coding paths and variables that might be subjected to changes.