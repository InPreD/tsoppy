# tsoppy

TSO500 v2 post processing cli

## Development

If you are contributing to the codebase of tsoppy, we would like for you to read this section to familiarize yourself with the structure and guidelines of the project.

### Recommended packages

The following packages are recommended to be used for the indicated purposes. We are open to discuss any of the following decisions depending on cost and benefit of the switch. If new purposes arise we will expand the list. Prior to suggesting new packages, please do some research on what people recommend and if the suggested package is actually maintained and fairly adopted by the community.

purpose | package
---|---
cli | [typer](https://typer.tiangolo.com/)
linting | [isort](https://pycqa.github.io/isort/), [ruff](https://docs.astral.sh/ruff/)
logging | [logging](https://docs.python.org/3/library/logging.html)
plotting | [plotnine](https://plotnine.org/)
read csv/tsv tables | [polars](https://docs.pola.rs/)
testing | [pytest](https://docs.pytest.org/en/stable/)
vcf | [CyVCF](https://github.com/arq5x/cyvcf/blob/master/README.rst)
