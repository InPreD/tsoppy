# Metric Plots Architecture

This document describes the internal architecture of the `metric_plots` module, its design decisions, and the responsibilities of each component. It is intended for developers extending or maintaining the metric plotting functionality.

## Overview

The `metric_plots` subpackage has two independent responsibilities:

1. **Extract and normalize QC metrics**
2. **Generate QC plots**

Separating these responsibilities allows metrics to be extracted once and reused for plotting, reporting, validation, or downstream analyses without reparsing workflow outputs.

The subpackage supports multiple TSO500 workflow implementations while exposing a single standardized data model to the plotting layer.

# High-level Architecture

```
                  Workflow Outputs
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
        ┌────────────────┴────────────────┐
        │                                 │
        ▼                                 ▼
master_metrics_table.tsv      joint_sequencing_QC_file.tsv
        │                                 │
        └────────────────┬────────────────┘
                         ▼
                 Polars DataFrames
                         │
                         ▼
                Python Plotting Module
                         │
                         ▼
                    QC Figures
```

Each layer has a single responsibility and communicates using standardized data structures.

---

# Design Principles

The module is built around the following principles:

- Single responsibility per component
- Automatic TSO500 workflow detection
- Standardized metric naming
- Workflow-independent plotting
- In-memory processing using Polars
- Minimal intermediate files
- Extensible architecture for future workflows

---

# Component Responsibilities

## MetricsOutputTsv

`MetricsOutputTsv` is responsible for parsing individual `MetricsOutput.tsv` files from TSO500 workflow output.

Responsibilities include:

- Reading `MetricsOutput.tsv`
- Parsing metric sections
- Extracting metadata
- Detecting workflow type and version
- Reading metric values
- Reading Lower Specification Limit (LSL) and Upper Specification Limit (USL) values

This class **does not perform plotting** and **does not rename metrics**. It simply represents one workflow output in a structured form.

---

## MetricPlots

`MetricPlots` converts parsed workflow outputs into standardized tables.

Responsibilities include:

- Iterating over workflow outputs
- Standardizing metric names
- Adding DNA/RNA prefixes
- Creating sample-level and threshold records
- Writing `master_metrics_table.tsv`and `joint_sequencing_QC_file.tsv`
- Producing Polars DataFrames for downstream plotting

This class forms the interface between parsing and visualization.

---

## Python Plotting Module

The plotting module is implemented entirely in Python.

Responsibilities include:

- Reading standardized Polars DataFrames
- Selecting runs for visualization
- Drawing QC plots
- Applying threshold lines
- Labelling plots
- Writing figures to file

The plotting code is intentionally isolated from the `MetricsOutput.tsv` parsing and only depends on the standardized tables.

---

# Data Model

The plotting module never needs to understand workflow-specific files.

Instead it operates on a normalized schema.

## Metadata

Each row contains metadata describing its origin.

Examples include:

- `RUN`
- `SAMPLE_ID`
- `WORKFLOW_TYPE`
- `WORKFLOW_VERSION`
- `RECORD_TYPE`

This allows multiple workflow versions to coexist in the same dataset.

---

## Metric Naming Convention

Metrics are standardized before plotting as follows:

| Prefix | Description |
|---------|-------------|
| *(none)* | Run-level metric |
| DNA_ | DNA-specific metric |
| RNA_ | RNA-specific metric |

This avoids ambiguity when DNA and RNA workflows expose metrics with identical names.

---

## Record Types

Each processed run generates multiple logical records.

```
DNA_SAMPLE

RNA_SAMPLE

LOWER_THRESHOLD

UPPER_THRESHOLD
```

Threshold records store the recommended specification limits reported by the workflow.

Keeping threshold values inside the standardized tables ensures plots always use the limits associated with the originating workflow version.

---

# Workflow Detection

Users are **not** required to specify:

- workflow type
- workflow version

Both are detected automatically from the workflow output.

This prevents accidental misclassification while allowing mixed workflow versions to be processed together.

---

# Polars DataFrames

Internally, all downstream processing uses Polars.

Reasons for choosing Polars include:

- efficient columnar execution
- low memory usage
- excellent performance on wide QC tables
- expressive dataframe API
- straightforward filtering of selected runs

The plotting module operates directly on Polars DataFrames.

No temporary plotting tables are written to disk.

---

# Plot Selection

Only a subset of runs may be plotted.

The module currently supports two selection strategies.

## Most recent runs

```
--plot-last-runs N
```

Selects the newest *N* runs from the master metrics table.

---

## Explicit run IDs

```
--plot-run-ids
```

Plots only the requested runs.

The two options are mutually exclusive.

---

# Why Threshold Rows?

Threshold values are stored as ordinary records rather than configuration files. Advantages include:

- version-specific thresholds remain attached to each workflow
- historical runs preserve historical limits
- plotting requires no external threshold configuration
- downstream analysis remains reproducible

---

# Supported Workflow Types

The current implementation supports:

- DRAGEN
- LocalApp

Both workflows are normalized into the same internal schema.

Adding support for a future workflow should require only parser extensions without modifying the plotting code.

---

# Extending the Parser

Adding a new metric typically requires:

1. Reading the metric from `MetricsOutput.tsv`
2. Mapping it to the standardized metric name
3. Adding it to the normalized output table

No plotting changes are required unless the metric should appear in a figure.

---

# Adding a New Plot

A new plot should:

1. Read the standardized Polars DataFrame.
2. Filter the required metric(s).
3. Use the associated threshold rows where applicable.
4. Produce the figure.

Plots should never parse workflow files directly.

---

# Error Handling

The parser is designed to tolerate differences between workflow versions.

Where possible:

- missing metrics are represented as null values
- workflow metadata is preserved
- unsupported metrics do not prevent extraction of supported metrics

This allows datasets from different workflow versions to coexist.

---

# Testing Strategy

Tests are divided into two logical groups.

## Parser tests

Verify that:

- workflow outputs are parsed correctly
- metadata is extracted correctly
- thresholds are detected correctly
- standardized tables are generated correctly

## Plotting tests

Verify that:

- selected runs are filtered correctly
- thresholds are applied correctly
- expected figures are generated

Keeping parser and plotting tests independent simplifies maintenance.

---

# Future Improvements

Potential future enhancements include:

- Interactive HTML QC dashboards
- Additional QC visualizations
- User-selectable plotting themes
- Statistical trend analysis across runs
- Automatic outlier detection
- Multi-instrument QC comparisons

The current architecture intentionally separates parsing from visualization so these features can be added without modifying the extraction pipeline.