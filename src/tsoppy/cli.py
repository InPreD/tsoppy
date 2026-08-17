"""This module contains the CLI commands for tsoppy."""

import importlib.metadata
import logging
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Annotated

import typer

from tsoppy.add_variant_summary.main import Summarise_variant_output
from tsoppy.metric_plots.main import MetricPlots
from tsoppy.metric_plots.plotting import Generate_qc_plots
from tsoppy.report_predispositions.main import generate_report

app = typer.Typer()

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s: %(message)s",
    datefmt="%Y/%m/%d %H:%M:%S",
)
logger = logging.getLogger(__name__)

app_version = importlib.metadata.version("tsoppy")


class WorkflowType(str, Enum):
    """Supported workflow types for metric plotting."""

    DRAGEN = "dragen"
    LOCALAPP = "localapp"


@dataclass
class CommonOptions:
    """Holds common input options."""

    nomenclature: Path
    config: Path


def _parse_run_ids(value: str) -> list[str]:
    """Split a comma-separated run ID string into stripped, non-empty IDs."""
    return [run_id.strip() for run_id in value.split(",") if run_id.strip()]


def _parse_run_id_file(value: str) -> list[str]:
    """Read run IDs from a file, one or more comma-separated IDs per line."""
    path = Path(value)

    if not path.is_file():
        raise typer.BadParameter(f"File '{value}' does not exist or is not a file.")

    return [
        run_id
        for line in path.read_text().splitlines()
        if line.strip()
        for run_id in _parse_run_ids(line)
    ]


@app.callback()
def main(
    ctx: typer.Context,
    config: Annotated[
        Path | None,
        typer.Option("--config", help="Path to tsoppy config file.", exists=True),
    ] = None,
    nomenclature: Annotated[
        Path | None,
        typer.Option("--nomenclature", help="Path to inpred nomenclature file."),
    ] = None,
):
    """
    This is the main entry point for the tsoppy CLI application. It sets up the
    context object with commonOptions that can be accessed by subcommands.
    """
    ctx.obj = CommonOptions(nomenclature=nomenclature, config=config)


@app.command()
def version():
    """Print the version of tsoppy."""
    print(f"tsoppy version {app_version}")


@app.command()
def metric_plots(
    input_glob: Annotated[
        str,
        typer.Option(
            help=(
                "Glob pattern matching workflow output directories "
                "whose final directory name is the sequencing run ID. "
                "Example:"
                " 'results/*/*' "
            ),
        ),
    ],
    inpred_nomenclature: Annotated[
        Path,
        typer.Option(
            exists=True,
            file_okay=True,
            dir_okay=False,
            readable=True,
            resolve_path=True,
            help="InPreD nomenclature YAML.",
        ),
    ],
    config_yaml: Annotated[
        Path,
        typer.Option(
            exists=True,
            file_okay=True,
            dir_okay=False,
            readable=True,
            resolve_path=True,
            help="Workflow configuration YAML.",
        ),
    ] = Path("config.yaml"),
    run_ids: Annotated[
        # Annotated as str, not list[str], because typer/click treats a
        # list[str]-annotated option as repeatable and would wrap this
        # parser's own list return value in another list. The annotation
        # only controls CLI parsing arity; the parser is what actually
        # decides the value this function receives.
        str | None,
        typer.Option(
            parser=_parse_run_ids,
            metavar="TEXT",
            help=(
                "Comma-separated run IDs to include in the generated master metrics "
                "table. If neither --run-ids nor --run-id-file is provided, all runs "
                "matched by --input-glob are included. Mutually exclusive with "
                "--run-id-file."
            ),
        ),
    ] = None,
    run_id_file: Annotated[
        # Stays annotated as str, not Path: Typer re-wraps a Path-annotated
        # option's parsed value in Path(...) again after the parser= runs,
        # which crashes since our parser returns a list[str]. _parse_run_id_file
        # does its own existence check instead of relying on typer.Option's
        # exists=True (which only applies to Path-annotated options).
        str | None,
        typer.Option(
            parser=_parse_run_id_file,
            metavar="FILE",
            help=(
                "Text file containing run IDs for generation of the master metrics "
                "table, one or more comma-separated per line. If neither "
                "--run-id-file nor --run-ids is provided, all runs matched by "
                "--input-glob are included. Mutually exclusive with --run-ids."
            ),
        ),
    ] = None,
    plot_run_ids: Annotated[
        # See the run_ids option above for why this stays annotated as str.
        str | None,
        typer.Option(
            parser=_parse_run_ids,
            metavar="TEXT",
            help=(
                "Comma-separated list of run IDs to include in plot. "
                "Mutually exclusive with --plot-run-id-file."
            ),
        ),
    ] = None,
    plot_run_id_file: Annotated[
        # See the run_id_file option above for why this stays annotated as str.
        str | None,
        typer.Option(
            parser=_parse_run_id_file,
            metavar="FILE",
            help=(
                "Text file containing run IDs to select for plotting, one or more "
                "comma-separated per line. Mutually exclusive with --plot-run-ids."
            ),
        ),
    ] = None,
    plot_last_runs: Annotated[
        int | None,
        typer.Option(
            min=1,
            help=(
                "Plot the most recent N runs for the selected workflow. "
                "Mutually exclusive with --plot-run-ids and "
                "--plot-run-id-file."
            ),
        ),
    ] = None,
    plot_workflow: Annotated[
        WorkflowType | None,
        typer.Option(
            help=(
                "Workflow whose runs will be plotted. If no plot run selector is "
                "provided, the last 10 runs for this workflow are plotted."
            ),
        ),
    ] = None,
    plot_output: Annotated[
        Path | None,
        typer.Option(
            help=(
                "Path to write the QC plot PDF to. Defaults to "
                "'<plot-workflow>_metric_plots.pdf' in the current directory."
            ),
        ),
    ] = None,
):
    """Create metrics tables and optionally generate QC plots."""
    logger.info("Creating metrics master table and joint QC.")

    # The mutual-exclusivity checks below all run here in the function body
    # rather than as per-option callbacks: Click invokes option callbacks in
    # command-line argument order, not declaration order, so a callback
    # checking ctx.params for another option could silently miss a conflict
    # depending on the order options are typed on the command line.
    if run_ids is not None and run_id_file is not None:
        message = "--run-id-file cannot be used with --run-ids."
        logger.error(message)
        raise typer.BadParameter(message)

    if plot_run_ids is not None and plot_run_id_file is not None:
        message = "--plot-run-id-file cannot be used with --plot-run-ids."
        logger.error(message)
        raise typer.BadParameter(message)

    plot_run_selection_given = plot_run_ids is not None or plot_run_id_file is not None

    prepare_plot_frames = (
        plot_workflow is not None
        or plot_last_runs is not None
        or plot_run_selection_given
    )

    # --plot-last-runs is a separate selection mode.
    if plot_last_runs is not None and plot_run_selection_given:
        message = (
            "--plot-last-runs cannot be combined with "
            "--plot-run-ids or --plot-run-id-file."
        )
        logger.error(message)
        raise typer.BadParameter(message)

    if prepare_plot_frames and plot_workflow is None:
        message = "--plot-workflow is required when plotting is requested."
        logger.error(message)
        raise typer.BadParameter(message)

    # run_ids/run_id_file are already parsed into a list by the options'
    # parser=, if given.
    resolved_run_ids = run_ids if run_ids is not None else run_id_file

    metric_plotter = MetricPlots(
        config_yaml=config_yaml,
        inpred_nomenclature=inpred_nomenclature,
        input_glob=input_glob,
        run_ids=resolved_run_ids,
    )

    master, joint_qc = metric_plotter.generate_metrics_tables()

    logger.info("Metrics master table and joint QC files created.")

    if prepare_plot_frames:
        # plot_run_ids/plot_run_id_file are already parsed into a list by the
        # options' parser=, if given.
        plotting_run_ids = (
            plot_run_ids if plot_run_ids is not None else plot_run_id_file
        )

        if plotting_run_ids is not None:
            # dict.fromkeys deduplicates while preserving order.
            plotting_run_ids = list(dict.fromkeys(plotting_run_ids))

        if plot_workflow is None:
            logger.error(
                "Internal error: --plot-workflow was not resolved before "
                "plot selection despite passing the earlier validation."
            )
            raise RuntimeError

        plot_frame, plot_joint_qc = metric_plotter.select_plot_data(
            master=master,
            joint_qc=joint_qc,
            workflow_type=plot_workflow.value,
            plot_last_runs=plot_last_runs,
            plot_run_ids=plotting_run_ids,
        )

        logger.info(
            f"Prepared {plot_frame.height} metric rows and "
            f"{plot_joint_qc.height} joint QC rows for {plot_workflow.value} plotting."
        )

        output_pdf = plot_output or Path(f"{plot_workflow.value}_metric_plots.pdf")

        Generate_qc_plots(
            metrics_table=plot_frame,
            joint_qc_table=plot_joint_qc,
            workflow=plot_workflow.value,
            output_pdf=output_pdf,
        )


@app.command()
def add_variant_summary(
    config_yaml: Annotated[
        Path, typer.Option("--config-yaml", "-c", help="Path to the config file.")
    ],
    inpred_nomenclature: Annotated[
        Path,
        typer.Option(
            "--inpred-nomenclature", "-i", help="Path to the inpred nomenclature file."
        ),
    ],
    results_directory: Annotated[
        Path,
        typer.Option(
            "--results-directory",
            "-r",
            help="path to the directory containing *_CombinedVariantOutput.tsv files.",
        ),
    ],
    output_file: Annotated[
        Path, typer.Option("--output-file", "-o", help="Path to the output file.")
    ],
):
    """Parse *_CombinedVariantOutput.tsv files into a variant summary table."""
    Summarise_variant_output(
        config_yaml, inpred_nomenclature, results_directory, output_file
    )



@app.command()
def report_predispositions(
    sample_id: Annotated[str, typer.Option(help="ID of the input sample.")],
    config_yaml: Annotated[str | Path, typer.Option(help="Config yaml file.")],
    inpred_nomenclature: Annotated[
        str | Path, typer.Option(help="InPreD nomenclature definition file.")
    ],
    root_path: Annotated[
        str | Path,
        typer.Option(
            help="Root path to the folder containing TSO500 analysis the sequencing run which contains the input sample."
        ),
    ],
    version_string: Annotated[
        str,
        typer.Option(
            help="String of the tools (and their versions) that were used in the data analysis."
        ),
    ],
    length_of_targeted_coding_regions: Annotated[
        float,
        typer.Option(
            help="Cummulative length of all coding regions targeted in the TSO500 experiment (in millions of bases)."
        ),
    ],
    tumor_purity: Annotated[
        float,
        typer.Option(
            help="Proportion of cancerous cells in the analysed sample, the value has to be in the range between 0 and 1, where 1 represents 100%."
        ),
    ],
    cancer_susceptibility_genes: Annotated[
        Path,
        typer.Option(
            help="File containing info about genes in which variants increase a person's risk of developing certain cancers."
        ),
    ],
    cancer_susceptibility_genes_column_list: Annotated[
        list[str],
        typer.Option(
            help="List of column names for the columns that should be reported in the output predisposition file."
        ),
    ],
    cancer_susceptibility_genes_gene_name_column: Annotated[
        str,
        typer.Option(
            help="Name of the column containing HUGO gene names (e.g. BRCA1, TP53)."
        ),
    ],
    output: Annotated[
        Path, typer.Option(help="File containing predisposition variants.")
    ],
):
    """
    This function reports variants called by small variant caller that are present
    in the cancer susceptibility genes.
    """

    generate_report(
        sample_id,
        config_yaml,
        inpred_nomenclature,
        root_path,
        version_string,
        length_of_targeted_coding_regions,
        tumor_purity,
        cancer_susceptibility_genes,
        cancer_susceptibility_genes_column_list,
        cancer_susceptibility_genes_gene_name_column,
        output,
    )
