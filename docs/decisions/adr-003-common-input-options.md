# Common input options

## Context and Problem Statement

We have input options that we would like to make available to most subcommands, e.g. `--config-yaml` and `--nomenclature-yaml`, without redefining them each time. Redefinition in each subcommand makes it easier to make mistakes and might introduce slight differences between subcommands.

### Option 1

Using typer context to propagate common input options to all relevant subcommands:

```python
import typer
from typing_extensions import Annotated
from dataclasses import dataclass

app = typer.Typer()

@dataclass
class CLIConfig:
    """Holds global configuration options."""
    env: str
    debug: bool

@app.callback()
def main(
    ctx: typer.Context,
    env: Annotated[str, typer.Option("--env", help="Target environment.")] = "production",
    debug: Annotated[bool, typer.Option("--debug", help="Enable debug output.")] = False
):
    ctx.obj = CLIConfig(env=env, debug=debug)

@app.command()
def deploy(ctx: typer.Context):
    # Access global values cleanly through the local context injection
    config = ctx.obj

    print(f"Deploying to environment: {config.env}")
    if config.debug:
        print("Debugging pipeline initialized...")
```

the subcommand could then be called with:

```bash
# Correct usage
$ python main.py --env staging deploy
# Incorrect usage (will throw an error)
$ python main.py deploy --env staging
```

The input options are defined once and will be applied to all subcommands that take the context. Calling the subcommand might seem a bit counterintuitive with flags being set both before and after the subcommand.

### Option 2

Define typer option globally and reuse the definition:

```python
import typer

app = typer.Typer()

# Define common options once to keep configurations DRY
VERBOSE_OPTION = typer.Option(False, "--verbose", "-v", help="Enable verbose logging.")
OUTPUT_OPTION = typer.Option("json", "--format", "-f", help="Output formatting style.")

@app.command()
def export_data(
    verbose: bool = VERBOSE_OPTION,
    output_format: str = OUTPUT_OPTION,
):
    typer.echo(f"Exporting... Format: {output_format}, Verbose: {verbose}")

@app.command()
def import_data(
    file_path: str,
    verbose: bool = VERBOSE_OPTION,
    output_format: str = OUTPUT_OPTION,
):
    typer.echo(f"Importing {file_path}... Format: {output_format}, Verbose: {verbose}")
```

the subcommand could then be called with:

```bash
$ python main.py import_data --file-path test.txt --verbose --output_format "yaml"
```

The input options are defined as global typer options but the correct option name must be added to the subcommand to make it consistent - more error prone. Calling the subcommand might feel more natural with all flags following the subcommand.

## Decision

## Consequences

### Positive

-

### Negative

-

### Neutral

-