# Metric plots

<<<<<<< Updated upstream
The `metric-plots` command extracts quality control (QC) metrics from TSO500 workflow outputs and produces standardized metrics tables together with QC plots.

The command supports both **DRAGEN** and **LocalApp** workflow outputs. Workflow type and workflow version are detected automatically from each workflow output, allowing multiple workflow versions to be processed simultaneously.
=======
The `metric-plots` subcommand collects quality-control (QC) metrics from TSO500 workflow outputs, standardizes them, writes consolidated metrics tables, and selects workflow-specific data for Python QC plotting.

The command supports **DRAGEN** and **LocalApp** workflow outputs. Workflow type and workflow version are detected from each workflow output.
>>>>>>> Stashed changes

For implementation details, see [`docs/references/metric_plots_architecture.md`](../references/metric_plots_architecture.md).

> **Current plotting status**
>
> The current implementation generates the standardized tables and prepares the in-memory `plot_frame` and `plot_joint_qc` DataFrames. The final call to the Python plotting function is still pending integration in `cli.py`.

## Quick start

Generate metrics tables for three runs:

```bash
tsoppy metric-plots \
    --input-glob 'results/*/*' \
    --run-ids RUN001,RUN002,RUN003
```

Generate the tables and select the last three DRAGEN runs for plotting:

```bash
tsoppy metric-plots \
    --input-glob 'results/*/*' \
    --run-ids RUN001,RUN002,RUN003 \
    --plot-workflow dragen \
    --plot-last-runs 3
```

View all available options with:

```bash
tsoppy metric-plots --help
```

## Input layout

`--input-glob` must match workflow output directories whose final directory name is the sequencing run ID.

Example:

```text
results/
├── dragen/
│   ├── RUN001/
│   ├── RUN002/
│   └── RUN003/
└── localapp/
    ├── RUN001/
    ├── RUN002/
    └── RUN003/
```

Use:

```bash
--input-glob 'results/*/*'
```

Quote the glob so that the shell passes it unchanged to `tsoppy`.

The glob must match workflow roots, not `MetricsOutput.tsv` files. A selected run is matched by exact directory name. For example, `RUN01` does not match `RUN010`.

A run may have both a DRAGEN and a LocalApp output. When both are valid and matched by the glob, both are processed.

## Selecting runs

Runs included in the generated metrics tables must be supplied using exactly one of:

- `--run-ids`
- `--run-id-file`

The two options are mutually exclusive.

### Comma-separated run IDs

```bash
tsoppy metric-plots \
    --input-glob 'results/*/*' \
    --run-ids RUN001,RUN002,RUN003
```

### Run ID file

```bash
tsoppy metric-plots \
    --input-glob 'results/*/*' \
    --run-id-file run_ids.txt
```

Example `run_ids.txt`:

```text
RUN001
RUN002
RUN003
```

Blank lines and lines beginning with `#` are ignored. A line may also contain comma-separated run IDs.

Duplicate IDs are removed while preserving their first occurrence.

## Run order and `RUN_INDEX`

Run order is significant.

The command does not parse dates from run identifiers. The last supplied run is treated as the latest/current run and receives:

```text
RUN_INDEX = 001
```

Earlier runs receive increasing values:

```text
RUN001 -> 003
RUN002 -> 002
RUN003 -> 001
```

## Configuration files

The default configuration files are:

```text
config.yaml
resources/nomenclature.yaml (to be changed/adjusted)
```

Alternative paths can be supplied with:

```bash
--config-yaml /path/to/config.yaml
--inpred-nomenclature /path/to/nomenclature.yaml
```

Both files must exist and be readable.

## Workflow detection

Workflow type and version are detected automatically and stored in:

```text
WORKFLOW_TYPE
WORKFLOW_VERSION
```

Users do not provide these values when generating the metrics tables.

`--plot-workflow` only selects which detected workflow type should be prepared for plotting. Supported values are: `dragen` and `localapp`

## Selecting runs for plotting

When plot selection is requested, `--plot-workflow` is required.

Choose one of:

- `--plot-last-runs`
- `--plot-run-ids`
- `--plot-run-id-file`

These three options are mutually exclusive.

### Last N workflow runs

```bash
tsoppy metric-plots \
    --input-glob 'results/*/*' \
    --run-id-file run_ids.txt \
    --plot-workflow dragen \
    --plot-last-runs 10
```

The command first filters to DRAGEN, then selects the last ten available DRAGEN runs in the preserved run order.

`--plot-last-runs` must be at least `1`. If fewer than `N` runs are available, all available runs are selected.

### Explicit plot run IDs

```bash
tsoppy metric-plots \
    --input-glob 'results/*/*' \
    --run-ids RUN001,RUN002,RUN003,RUN004 \
    --plot-workflow localapp \
    --plot-run-ids RUN002,RUN004
```

Requested runs that are unavailable for the selected workflow are skipped with a warning.

### Plot run ID file

```bash
tsoppy metric-plots \
    --input-glob 'results/*/*' \
    --run-id-file run_ids.txt \
    --plot-workflow localapp \
    --plot-run-id-file plot_run_ids.txt
```

Example:

```text
RUN002
RUN004
```

Blank lines and comment lines are ignored.

## CLI options

| Option | Requirement | Purpose |
|---|---|---|
| `--input-glob` | Required | Match workflow root directories |
| `--config-yaml` | Optional | Workflow configuration file |
| `--inpred-nomenclature` | Optional | InPreD nomenclature file |
| `--run-ids` | Exactly one master selector | Comma-separated run IDs |
| `--run-id-file` | Exactly one master selector | File containing run IDs |
| `--plot-workflow` | Required with plot selection | Select `dragen` or `localapp` |
| `--plot-last-runs` | Optional plot selector | Select the last `N` workflow runs |
| `--plot-run-ids` | Optional plot selector | Select explicit plot run IDs |
| `--plot-run-id-file` | Optional plot selector | Select plot runs from a file |

## Generated outputs

The command writes two canonical files to the current working directory:

```text
master_metrics_table.tsv
joint_sequencing_QC_file.tsv
```

There is no `--workdir` or output-directory option. Workflow managers such as Nextflow are expected to manage the process working directory and publish outputs afterwards.

### `master_metrics_table.tsv`

The master table contains standardized sample-level, threshold, and run-level metrics.

Metadata columns are:

| Column | Description |
|---|---|
| `RUN_INDEX` | Run-order index; `001` is the last supplied run |
| `SAMPLE_ID` | Sample or threshold identifier |
| `RUN` | Sequencing run ID |
| `WORKFLOW_TYPE` | Detected workflow type |
| `WORKFLOW_VERSION` | Detected workflow version |
| `RECORD_TYPE` | Logical row type |

Metric columns follow this naming convention:

| Prefix | Meaning |
|---|---|
| none | Run-level or workflow-wide metric |
| `DNA_` | DNA-specific metric |
| `RNA_` | RNA-specific metric |

Finalized missing values are represented as `NA`.

Record types include:

```text
DNA_SAMPLE
RNA_SAMPLE
SAMPLE
LOWER_THRESHOLD
UPPER_THRESHOLD
```

Sample type is determined primarily from the presence of DNA- or RNA-specific metric values. Explicit `DNA_` or `RNA_` sample ID prefixes are used only as a fallback.

Threshold rows keep lower and upper specification guidelines associated with the workflow type and version that produced them.

### `joint_sequencing_QC_file.tsv`

The joint QC table contains one row per run and workflow.

Columns are:

```text
RUN_INDEX
RUN_ID
WORKFLOW_TYPE
WORKFLOW_VERSION
PCT_PF_READS
PCT_Q30_R1
PCT_Q30_R2
CLUSTER_DENSITY
ESTIMATED_YIELD
CLUSTERS_PASSING_FILTER
```

Metrics unavailable for a workflow are represented as `NA`.

## In-memory plotting data

`select_plot_data()` returns:

```text
plot_frame
plot_joint_qc
```

`plot_frame` contains the selected rows from the master table.

`plot_joint_qc` contains matching sequencing-QC rows.

These DataFrames are intended to be passed directly to the plotting implementation.

## Complete examples

Generate tables only:

```bash
tsoppy metric-plots \
    --input-glob '/data/tso500/*/*' \
    --run-id-file run_ids.txt
```

Select the last eight DRAGEN runs:

```bash
tsoppy metric-plots \
    --input-glob '/data/tso500/*/*' \
    --run-id-file run_ids.txt \
    --plot-workflow dragen \
    --plot-last-runs 8
```

Select explicit LocalApp runs:

```bash
tsoppy metric-plots \
    --input-glob '/data/tso500/*/*' \
    --run-ids RUN001,RUN002,RUN003,RUN004 \
    --plot-workflow localapp \
    --plot-run-ids RUN002,RUN004
```

Use custom configuration files:

```bash
tsoppy metric-plots \
    --input-glob '/data/tso500/*/*' \
    --config-yaml /data/config.yaml \
    --inpred-nomenclature /data/nomenclature.yaml \
    --run-id-file run_ids.txt
```

## Troubleshooting

### The input glob matches nothing

Check that it points to workflow directories, not directly to metrics files:

```bash
--input-glob 'results/*/*'
```

The final directory name must exactly equal the requested run ID.

### A requested plot run is skipped

The run may be part of the master run list but unavailable for the selected workflow. The command logs a warning and continues with available requested runs.

### Plot workflow is missing

When using a plot selector, add:

```bash
--plot-workflow dragen
```

or:

```bash
--plot-workflow localapp
```

For internal processing and extension guidance, see [`docs/references/metric_plots_architecture.md`](../references/metric_plots_architecture.md).
