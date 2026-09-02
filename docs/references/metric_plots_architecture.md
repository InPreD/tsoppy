# Metric plots architecture

This document describes the high-level architecture of the `metric_plots` functionality.

For command-line usage, options, and examples, see [`docs/guides/metric_plots.md`](../guides/metric_plots.md).

## Design

The implementation separates workflow discovery, metrics parsing, normalization, and plotting.

**Parse workflow-specific data once, normalize it once, and let downstream consumers operate on the standardized representation.**

```mermaid
%%{init: {'flowchart': {'defaultRenderer': 'diagre-wrapper'}} }%%
flowchart TD
    CLI(CLI) --> WorkflowOutput(WorkflowOutput)
    WorkflowOutput --> MetricsOutputTsv(MetricsOutputTsv)
    MetricsOutputTsv --> MetricPlots(MetricPlots)
    MetricPlots --> Tables(Standardized metrics tables)
    Tables --> PlotSelection(Plot-data selection)
    PlotSelection --> Plotting(Plotting)
```

## Components

### CLI

The CLI validates user input, creates `MetricPlots`, generates the metrics tables, and optionally requests plotting.

Metric transformation logic is kept outside the CLI.

### `WorkflowOutput`

`WorkflowOutput` represents a workflow output directory and provides the workflow-aware context needed to locate downstream files.

### `MetricsOutputTsv`

`MetricsOutputTsv` parses individual `MetricsOutput.tsv` files and exposes their metric sections together with workflow type and version.

### `MetricPlots`

`MetricPlots` normalizes metrics from different TSO500 workflows into a common representation.

It is responsible for combining selected runs, assigning run metadata, creating the master metrics table and joint sequencing QC table, and selecting data for plotting.

### Plotting

The plotting layer consumes standardized Polars DataFrames.

It does not discover workflow outputs or parse `MetricsOutput.tsv` files.

Plot definitions are kept separate from the normalized metric data so the same plotting code can be used for supported workflow types.

## Standardized data

Metrics from DRAGEN and LocalApp are normalized into a common schema before plotting.

The master table contains metadata including:

- `RUN_INDEX`
- `SAMPLE_ID`
- `RUN`
- `WORKFLOW_TYPE`
- `WORKFLOW_VERSION`
- `RECORD_TYPE`

Workflow-specific DNA and RNA metrics use `DNA_` and `RNA_` prefixes.

Sample type is taken from `Sample_Type` in the workflow sample sheet. The sample sheet is the authoritative source because it determines how the upstream analysis is performed.

Lower and upper specification limits are represented as threshold records in the standardized data.

## Outputs

The command writes two canonical metrics files:

- `master_metrics_table.tsv`
- `joint_sequencing_QC_file.tsv`

When plotting is requested, the standardized in-memory tables are passed directly to the plotting layer. Temporary plotting tables are not written to disk.

## Run handling

Runs may be selected explicitly or discovered from the input glob.

Run ordering is resolved before metrics from individual workflow outputs are combined. `RUN_INDEX` records that resolved order and does not depend on filesystem glob ordering.

Workflow-specific plotting operates on the standardized tables rather than rereading workflow outputs.

Detailed run-selection behavior belongs in the [metric plots user guide](../guides/metric_plots.md).

## Validation

Validation is performed at the layer responsible for the corresponding input:

- Typer validates command-line values and file options.
- `MetricPlots` validates run selection and normalized data.
- workflow parser classes validate workflow-specific files.

Recoverable conditions are logged while invalid or ambiguous data that would make the output unreliable causes processing to fail.

## Testing

Metric processing, CLI behavior, and plotting are tested independently.

Fixtures should remain small and generic while covering supported workflows and important edge cases.

## Extension

New workflows should be integrated at the workflow/parsing layer and normalized into the existing standardized representation.

New metrics should enter the standardized tables before being consumed by downstream plotting.

New plots should consume standardized plotting data and should not parse raw workflow files or recreate upstream metadata.