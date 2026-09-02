# Metric plots

The `metric-plots` command extracts quality control (QC) metrics from TSO500 workflow outputs and produces standardized metrics tables together with QC plots.

The command supports both **DRAGEN** and **LocalApp** workflow outputs. Workflow type and workflow version are detected automatically from each workflow output, allowing multiple workflow versions to be processed simultaneously.

For implementation details, see [`docs/references/metric_plots_architecture.md`](../references/metric_plots_architecture.md).

## Quick start

Generate metrics tables for three runs:

```bash
tsoppy metric-plots \
    --input-glob 'results/*/*' \
    --inpred-nomenclature /path/to/nomenclature.yaml \
    --run-ids RUN001,RUN002,RUN003
```

Generate tables for five runs and plot the last three DRAGEN runs among them:

```bash
tsoppy metric-plots \
    --input-glob 'results/*/*' \
    --inpred-nomenclature /path/to/nomenclature.yaml \
    --run-ids RUN001,RUN002,RUN003,RUN004,RUN005 \
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

The glob must match workflow root directories, not `MetricsOutput.tsv` files. A specific run can be selected by pointing the glob at its workflow root directory, e.g. `results/localapp/RUN01`. Matching is by exact directory name, so `RUN01` does not match `RUN010`.

A run may have both a DRAGEN and a LocalApp output. When both are valid and matched by the glob, both are processed.

## Selecting runs

Runs included in the generated metrics tables can be supplied using at most one of:

- `--run-ids`
- `--run-id-file`

The two options are mutually exclusive. If neither is provided, all runs matched by `--input-glob` are included.

### Comma-separated run IDs

```bash
tsoppy metric-plots \
    --input-glob 'results/*/*' \
    --inpred-nomenclature /path/to/nomenclature.yaml \
    --run-ids RUN001,RUN002,RUN003
```

### Run ID file

```bash
tsoppy metric-plots \
    --input-glob 'results/*/*' \
    --inpred-nomenclature /path/to/nomenclature.yaml \
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

Run order is determined by the resolved run ID order, not by the order in which directories are returned by `--input-glob`.

For `--run-ids`, the comma-separated order is preserved. For `--run-id-file`, run IDs are processed in file order. If neither option is supplied, all run IDs matched by `--input-glob` are used and sorted lexicographically.

Duplicate run IDs are removed while preserving their first occurrence.

The command does not parse dates from run identifiers. The last run in the resolved order is treated as the latest/current run and receives:

```text
RUN_INDEX = 001
```

Earlier runs receive increasing values:

```text
RUN001 -> 003
RUN002 -> 002
RUN003 -> 001
```

## Workflow selection

`--plot-workflow` selects which workflow type should be prepared for plotting. Supported values are: `dragen` and `localapp`.

## Selecting runs for plotting

`--plot-workflow` is required to create plots.

If `--plot-workflow` is provided without a plot run selector, the last 10 available runs for that workflow are selected by default.

To override the default selection, use one of:

- `--plot-last-runs`
- `--plot-run-ids`
- `--plot-run-id-file`

These three options are mutually exclusive.

### Last N workflow runs

```bash
tsoppy metric-plots \
    --input-glob 'results/*/*' \
    --inpred-nomenclature /path/to/nomenclature.yaml \
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
    --inpred-nomenclature /path/to/nomenclature.yaml \
    --run-ids RUN001,RUN002,RUN003,RUN004 \
    --plot-workflow localapp \
    --plot-run-ids RUN002,RUN004
```

Requested runs that are unavailable for the selected workflow are skipped with a warning.

### Plot run ID file

```bash
tsoppy metric-plots \
    --input-glob 'results/*/*' \
    --inpred-nomenclature /path/to/nomenclature.yaml \
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
| `--input-glob` | Required | Glob pattern matching workflow output directories whose final directory name is the sequencing run ID |
| `--inpred-nomenclature` | Required | InPreD nomenclature YAML |
| `--config-yaml` | Optional; default: `config.yaml` | Workflow configuration YAML |
| `--run-ids` | Optional; at most one of `--run-ids` or `--run-id-file`; defaults to all runs matched by `--input-glob` | Comma-separated run IDs to include in the generated master metrics table |
| `--run-id-file` | Optional; at most one of `--run-ids` or `--run-id-file`; defaults to all runs matched by `--input-glob` | Text file containing run IDs for generation of the master metrics table |
| `--plot-run-ids` | Optional plot selector | Comma-separated run IDs to include in the plot |
| `--plot-run-id-file` | Optional plot selector | Text file containing run IDs to select for plotting |
| `--plot-last-runs` | Optional plot selector; integer ≥ 1 | Plot the most recent `N` runs for the selected workflow |
| `--plot-workflow` | Required when plotting is requested | Workflow to plot: `dragen` or `localapp` |
| `--help` | Optional | Show the command help and exit |

The plot selectors `--plot-run-ids`, `--plot-run-id-file`, and `--plot-last-runs` are mutually exclusive.

## Generated outputs

The command always writes two canonical files to the current working directory:

```text
master_metrics_table.tsv
joint_sequencing_QC_file.tsv
```

When plotting is requested with --plot-workflow, it also writes a workflow-specific PDF:

dragen_metric_plots.pdf or localapp_metric_plots.pdf

There is no --workdir or output-directory option. Workflow managers such as Nextflow are expected to manage the process working directory and publish outputs afterwards.
### `master_metrics_table.tsv`

The master table contains standardized sample-level, threshold, and run-level metrics.

Metadata columns are:

| Column | Description |
|---|---|
| `RUN_INDEX` | Run-order index; `001` is the last run in the resolved run order |
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

Sample type is determined from `Sample_Type` in the workflow SampleSheet. Metrics samples are matched to the SampleSheet using `Pair_ID` when available, otherwise `Sample_ID`. `DNA` and `RNA` values are assigned `DNA_SAMPLE` and `RNA_SAMPLE`, respectively. Samples without an unambiguous supported SampleSheet classification are assigned `SAMPLE`.

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

These DataFrames are passed directly to `Generate_qc_plots()`, which renders them to the output PDF.

## Complete examples

Generate tables only:

```bash
tsoppy metric-plots \
    --input-glob '/data/tso500/*/*' \
    --inpred-nomenclature /path/to/nomenclature.yaml \
    --run-id-file run_ids.txt
```

Select the last eight DRAGEN runs:

```bash
tsoppy metric-plots \
    --input-glob '/data/tso500/*/*' \
    --inpred-nomenclature /path/to/nomenclature.yaml \
    --run-id-file run_ids.txt \
    --plot-workflow dragen \
    --plot-last-runs 8
```

Select explicit LocalApp runs:

```bash
tsoppy metric-plots \
    --input-glob '/data/tso500/*/*' \
    --inpred-nomenclature /path/to/nomenclature.yaml \
    --run-ids RUN001,RUN002,RUN003,RUN004 \
    --plot-workflow localapp \
    --plot-run-ids RUN002,RUN004
```

Use custom configuration files:

```bash
tsoppy metric-plots \
    --input-glob '/data/tso500/*/*' \
    --inpred-nomenclature /path/to/nomenclature.yaml \
    --config-yaml /data/config.yaml \
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
