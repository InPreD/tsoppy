"""
This module contains the CLI commands for tsoppy.
"""

import importlib.metadata
import logging
from typing import Annotated

import typer

# Set up logging for the CLI. The logging level is set to INFO, and the log messages will include the timestamp, log level, and message.
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(levelname)s: %(message)s',
                    datefmt='%Y/%m/%d %H:%M:%S')
logger = logging.getLogger(__name__)

app = typer.Typer()
app_version = importlib.metadata.version("tsoppy")


@app.command()
def version():
    """
    Prints the version of tsoppy.
    """
    print(f"tsoppy version {app_version}")


@app.command()
def placeholder(
    user_name: Annotated[str, typer.Option("--name", "-n")],
    user_id: Annotated[str, typer.Option("--id", "-i")],
    verbose: Annotated[bool, typer.Option("--verbose", "-v")] = False,
):
    """
    This is the helptext for the placeholder command that demonstrates how to
    use Typer for CLI applications.
    """
    if verbose:
        print(f"{user_name} has the following id: {user_id}")
    else:
        print(f"{user_name}: {user_id}")
