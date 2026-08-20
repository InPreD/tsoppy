# Metric plots

The `metric-plots` command extracts quality control (QC) metrics from TSO500 workflow outputs and produces standardized metrics tables together with QC plots.

The command supports both **DRAGEN** and **LocalApp** workflow outputs. Workflow type and workflow version are detected automatically from each workflow output, allowing multiple workflow versions to be processed simultaneously.

---

# Overview

The workflow consists of two independent stages:

1. **Metrics extraction**
    - Reads workflow outputs
    - Parses MetricsOutput.tsv
    - Standardizes metric names
    - Writes consolidated metrics tables

2. **Plot generation**
    - Reads the standardized metrics tables directly in memory
    - Generates QC plots using the Python plotting module
    - No temporary plotting TSV files are created

Separating these stages allows the same metrics tables to be reused for additional analyses without rerunning the parser.

---

# Input

The command expects an input directory containing workflow result folders.

Supported workflow outputs include

- DRAGEN
- LocalApp

The input directory may contain

- only DRAGEN runs
- only LocalApp runs
- a mixture of both

The workflow type and workflow version are detected automatically.

---

# Selecting runs

Run IDs can be supplied in two ways.

## Command line

```bash
tsoppy metric-plots \
    --input-directory /path/to/results \
    --run-ids RUN001,RUN002,RUN003
```

---

## Text file

One run ID per line

```text
RUN001
RUN002
RUN003
```

Run

```bash
tsoppy metric-plots \
    --input-directory /path/to/results \
    --run-id-file run_ids.txt
```

---

## Combine both

Both options may be supplied simultaneously.

Duplicate run IDs are removed automatically.

```bash
tsoppy metric-plots \
    --input-directory results \
    --run-id-file run_ids.txt \
    --run-ids RUN005,RUN006
```

---

# Output directory

Output is written to the current working directory by default.

A different directory may be specified using

```bash
--workdir
```

Example

```bash
tsoppy metric-plots \
    --input-directory results \
    --run-id-file run_ids.txt \
    --workdir output
```

---

# Output files

The command generates

```text
master_metrics_table.tsv

joint_sequencing_QC_file.tsv
```

These files become the canonical input for downstream QC plotting.

---

# master_metrics_table.tsv

This table contains all sample-level QC metrics collected from every processed workflow.

Metadata columns include

```text
RUN
SAMPLE_ID
WORKFLOW_TYPE
WORKFLOW_VERSION
RECORD_TYPE
```

Metric columns are standardized across workflow versions.

Naming convention

| Prefix | Meaning |
|---------|---------|
| none | Run-level metric |
| DNA_ | DNA metric |
| RNA_ | RNA metric |

Each processed workflow contributes multiple record types

```text
DNA_SAMPLE

RNA_SAMPLE

LOWER_THRESHOLD

UPPER_THRESHOLD
```

Threshold rows contain the LSL and USL values reported by the originating workflow.

Since workflow type and version are preserved, metrics originating from different workflow versions remain distinguishable.

> Metrics having identical names may represent different calculations in different workflow versions. Comparisons across workflow versions should therefore be interpreted carefully.

---

# joint_sequencing_QC_file.tsv

Contains sequencing-level metrics for each processed run.

Examples include

```text
PCT_PF_READS

PCT_Q30_R1

PCT_Q30_R2

CLUSTER_DENSITY

ESTIMATED_YIELD

CLUSTERS_PASSING_FILTER
```

This file is intended for sequencing QC reporting.

---

# Plot generation

QC plots are generated using the new Python plotting implementation.

The plotting module operates directly on Polars DataFrames produced during metrics extraction.

No intermediate plotting TSV files are written or reread.

This significantly reduces I/O and keeps the plotting pipeline entirely in memory.

---

# Plot selection

The plotting module supports two methods for selecting runs.

## Plot the latest runs

```bash
tsoppy metric-plots \
    --input-directory results \
    --run-id-file run_ids.txt \
    --plot-last-runs 8
```

The newest N runs are selected according to their ordering in the generated master metrics table.

All associated records are retained

- DNA_SAMPLE
- RNA_SAMPLE
- LOWER_THRESHOLD
- UPPER_THRESHOLD

---

## Plot selected runs

```bash
tsoppy metric-plots \
    --input-directory results \
    --plot-run-ids RUN001,RUN004,RUN010
```

The options

```text
--plot-last-runs

--plot-run-ids
```

are mutually exclusive.

---

# Internal workflow

```
Workflow outputs
        │
        ▼
 MetricsOutput.tsv
        │
        ▼
 MetricsOutputTsv
        │
        ▼
 MetricPlots
        │
        ├──────────────┐
        │              │
        ▼              ▼
master_metrics     joint sequencing QC
table              table
        │
        ▼
 Polars DataFrame
        │
        ▼
 Python plotting module
        │
        ▼
 QC figures
```

Responsibilities are intentionally separated.

**MetricsOutputTsv**

- Parses workflow MetricsOutput.tsv files
- Detects workflow type
- Detects workflow version

**MetricPlots**

- Converts parsed workflow metrics into standardized tables
- Produces master metrics table
- Produces sequencing QC table
- Creates filtered Polars DataFrames for plotting

**Python plotting module**

- Reads Polars DataFrames
- Generates QC plots
- Uses workflow metadata and threshold rows directly from the standardized tables

---

# Examples

Generate metrics only

```bash
tsoppy metric-plots \
    --input-directory results \
    --run-id-file run_ids.txt
```

Generate metrics and plot the latest ten runs

```bash
tsoppy metric-plots \
    --input-directory results \
    --run-id-file run_ids.txt \
    --plot-last-runs 10
```

Generate metrics and plot specific runs

```bash
tsoppy metric-plots \
    --input-directory results \
    --plot-run-ids RUN001,RUN002,RUN003
```

---

# Help

```bash
tsoppy metric-plots --help
```