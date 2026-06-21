#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import re
import subprocess
import sys
from collections import OrderedDict
from pathlib import Path
from typing import Iterable

import polars as pl


METRIC_HEADER_FIRST_COLUMN = "Metric (UOM)"
GUIDELINE_COLUMNS = {"LSL Guideline", "USL Guideline"}
BLANK_VALUES = {"", "NA", "N/A", "nan", "None"}
NON_SAMPLE_IDS = {"Value", "LSL Guideline", "USL Guideline", "Metric (UOM)"}
LOCALAPP_CANONICAL_COLUMNS = [
    "RUN_ID",
    "SAMPLE_ID",
    "PCT_PF_READS",
    "PCT_Q30_R1",
    "PCT_Q30_R2",
    "COMPLETED_ALL_STEPS",
    "FAILED_STEPS",
    "STEPS_NOT_EXECUTED",
    "CONTAMINATION_SCORE",
    "CONTAMINATION_P_VALUE",
    "MEDIAN_INSERT_SIZE",
    "MEDIAN_EXON_COVERAGE",
    "PCT_EXON_50X",
    "USABLE_MSI_SITES",
    "COVERAGE_MAD",
    "MEDIAN_BIN_COUNT_CNV_TARGET",
    "TOTAL_PF_READS",
    "MEAN_FAMILY_SIZE",
    "MEDIAN_TARGET_COVERAGE",
    "PCT_CHIMERIC_READS",
    "PCT_EXON_100X",
    "PCT_READ_ENRICHMENT",
    "PCT_USABLE_UMI_READS",
    "MEAN_TARGET_COVERAGE",
    "PCT_ALIGNED_READS",
    "PCT_CONTAMINATION_EST",
    "PCT_PF_UQ_READS",
    "PCT_TARGET_0_4X_MEAN",
    "PCT_TARGET_100X",
    "PCT_TARGET_250X",
    "MEDIAN_CV_GENE_500X",
    "TOTAL_ON_TARGET_READS",
    "PCT_ON_TARGET_READS",
    "SCALED_MEDIAN_GENE_COVERAGE",
]


def clean_cell(value: object) -> str:
    if value is None:
        return ""
    return str(value).strip()


def normalize_metric_name(metric_name: str) -> str:
    metric_name = clean_cell(metric_name)
    metric_name = re.sub(r"\s*\([^)]*\)\s*$", "", metric_name)
    metric_name = re.sub(r"%", "PCT", metric_name)
    metric_name = re.sub(r"[^0-9A-Za-z]+", "_", metric_name)
    metric_name = re.sub(r"_+", "_", metric_name).strip("_")
    return metric_name


def normalize_section_name(section: str) -> str:
    section = clean_cell(section).strip("[]")
    section = re.sub(r"[^0-9A-Za-z]+", "_", section)
    section = re.sub(r"_+", "_", section).strip("_")
    return section.upper()


def read_tsv_rows(path: Path) -> list[list[str]]:
    with path.open("r", newline="") as handle:
        return [[clean_cell(cell) for cell in row] for row in csv.reader(handle, delimiter="\t")]


def is_blank_row(row: list[str]) -> bool:
    return all(clean_cell(cell) == "" for cell in row)


def is_section_row(row: list[str]) -> bool:
    first = row[0] if row else ""
    return first.startswith("[") and first.endswith("]")


def pad_row(row: list[str], size: int) -> list[str]:
    if len(row) >= size:
        return row
    return row + [""] * (size - len(row))


def has_non_blank(values: Iterable[str]) -> bool:
    return any(clean_cell(value) != "" for value in values)


def is_missing(value: str) -> bool:
    return clean_cell(value) in BLANK_VALUES


def unique_metric_name(metric: str, section: str, existing: Iterable[str]) -> str:
    existing_set = set(existing)
    candidate = f"{normalize_section_name(section)}__{metric}" if section else metric
    if candidate not in existing_set:
        return candidate
    counter = 2
    while f"{candidate}_{counter}" in existing_set:
        counter += 1
    return f"{candidate}_{counter}"


def add_metric_value(sample_metrics, sample_id, metric, value, section, disambiguate_duplicates):
    sample_id = clean_cell(sample_id)
    metric = normalize_metric_name(metric)
    value = clean_cell(value)

    if sample_id == "" or sample_id in NON_SAMPLE_IDS or metric == "":
        return

    sample_metrics.setdefault(sample_id, OrderedDict())
    metrics = sample_metrics[sample_id]

    if metric not in metrics:
        metrics[metric] = value
        return

    existing = metrics[metric]
    if is_missing(existing) and not is_missing(value):
        metrics[metric] = value
        return

    if is_missing(value) or value == existing:
        return

    if disambiguate_duplicates:
        metrics[unique_metric_name(metric, section, metrics.keys())] = value


def add_run_metric(run_qc_metrics, metric, value, section, disambiguate_duplicates):
    metric = normalize_metric_name(metric)
    value = clean_cell(value)

    if metric == "":
        return

    if metric not in run_qc_metrics:
        run_qc_metrics[metric] = value
        return

    existing = run_qc_metrics[metric]
    if is_missing(existing) and not is_missing(value):
        run_qc_metrics[metric] = value
        return

    if is_missing(value) or value == existing:
        return

    if disambiguate_duplicates:
        run_qc_metrics[unique_metric_name(metric, section, run_qc_metrics.keys())] = value


def parse_metric_table(rows, start_index, current_section, sample_metrics, run_qc_metrics, disambiguate_duplicates, known_sample_ids):
    header = rows[start_index]
    value_col_idx = None
    sample_columns = []

    for col_idx, col_name in enumerate(header[1:], start=1):
        if col_name == "Value" and value_col_idx is None:
            value_col_idx = col_idx
        elif col_name not in GUIDELINE_COLUMNS and col_name not in NON_SAMPLE_IDS and col_name != "":
            sample_columns.append((col_idx, col_name))

    header_sample_names = [sample_id for _, sample_id in sample_columns]
    has_duplicate_sample_names = len(header_sample_names) != len(set(header_sample_names))
    data_column_count = max(len(header) - 3, 0)
    if known_sample_ids and has_duplicate_sample_names and len(known_sample_ids) <= data_column_count:
        sample_columns = [(idx, sample_id) for idx, sample_id in enumerate(known_sample_ids, start=3)]

    row_index = start_index + 1
    while row_index < len(rows):
        row = rows[row_index]
        if is_blank_row(row) or is_section_row(row):
            break

        row = pad_row(row, len(header))
        metric = row[0]

        if value_col_idx is not None and not sample_columns:
            add_run_metric(run_qc_metrics, metric, row[value_col_idx], current_section, disambiguate_duplicates)
        else:
            for col_idx, sample_id in sample_columns:
                value = row[col_idx] if col_idx < len(row) else ""
                add_metric_value(sample_metrics, sample_id, metric, value, current_section, disambiguate_duplicates)

        row_index += 1

    return row_index, [sample_id for _, sample_id in sample_columns]


def parse_analysis_status(rows, start_index, current_section, sample_metrics, disambiguate_duplicates):
    header = rows[start_index]
    first_data_index = start_index + 1
    if first_data_index >= len(rows):
        return first_data_index, []

    if header[0] == "":
        sample_columns = [(idx, sample_id) for idx, sample_id in enumerate(header) if sample_id and sample_id not in NON_SAMPLE_IDS]
    else:
        sample_columns = [(idx, sample_id) for idx, sample_id in enumerate(header, start=1) if sample_id and sample_id not in NON_SAMPLE_IDS]

    row_index = first_data_index
    while row_index < len(rows):
        row = rows[row_index]
        if is_blank_row(row) or is_section_row(row):
            break

        expected_size = max([idx for idx, _ in sample_columns], default=0) + 1
        row = pad_row(row, expected_size)
        metric = row[0]
        for col_idx, sample_id in sample_columns:
            value = row[col_idx] if col_idx < len(row) else ""
            add_metric_value(sample_metrics, sample_id, metric, value, current_section, disambiguate_duplicates)

        row_index += 1

    return row_index, [sample_id for _, sample_id in sample_columns]


def parse_metrics_output(metrics_path: Path, run_label: str, disambiguate_duplicates: bool = True) -> pl.DataFrame:
    rows = read_tsv_rows(metrics_path)
    sample_metrics = OrderedDict()
    run_qc_metrics = OrderedDict()
    current_section = ""
    row_index = 0
    known_sample_ids = []

    while row_index < len(rows):
        row = rows[row_index]

        if is_section_row(row):
            current_section = row[0].strip("[]")
            row_index += 1
            continue

        if current_section == "Analysis Status" and row and has_non_blank(row):
            row_index, analysis_sample_ids = parse_analysis_status(rows, row_index, current_section, sample_metrics, disambiguate_duplicates)
            if analysis_sample_ids:
                known_sample_ids = analysis_sample_ids
            continue

        if is_blank_row(row):
            row_index += 1
            continue

        if row and row[0] == METRIC_HEADER_FIRST_COLUMN:
            row_index, table_sample_ids = parse_metric_table(rows, row_index, current_section, sample_metrics, run_qc_metrics, disambiguate_duplicates, known_sample_ids)
            if not known_sample_ids and table_sample_ids:
                known_sample_ids = table_sample_ids
            continue

        row_index += 1

    if not sample_metrics:
        raise ValueError(f"No sample-wise metrics were parsed from {metrics_path}")

    records = []
    for sample_id, metrics in sample_metrics.items():
        if sample_id in NON_SAMPLE_IDS:
            continue
        record = OrderedDict()
        record["RUN_ID"] = run_label
        record["SAMPLE_ID"] = sample_id
        for metric, value in run_qc_metrics.items():
            record[metric] = value
        for metric, value in metrics.items():
            record[metric] = value
        records.append(record)

    return pl.DataFrame(records, infer_schema_length=None)


def make_unique_columns(columns: Iterable[str]) -> list[str]:
    counts = {}
    result = []
    for col in columns:
        counts[col] = counts.get(col, 0) + 1
        if counts[col] == 1:
            result.append(col)
        else:
            result.append(f"{col}_{counts[col] - 1}")
    return result


def run_id_from_filename(path: Path) -> str:
    name = path.name
    name = re.sub(r"_MetricsOutput_(Localapp|LocalApp|Dragen|DRAGEN)\.tsv$", "", name)
    name = re.sub(r"_MetricsOutput\.tsv$", "", name)
    return name


def read_run_ids(path: Path) -> list[str]:
    if not path.is_file():
        raise FileNotFoundError(f"RUN_IDs file not found: {path}")
    if path.suffix.lower() == ".csv":
        with path.open("r", newline="") as handle:
            rows = [[clean_cell(cell) for cell in row] for row in csv.reader(handle)]
    elif path.suffix.lower() == ".tsv":
        rows = read_tsv_rows(path)
    else:
        rows = [[line.strip()] for line in path.read_text().splitlines()]

    values = []
    for row in rows:
        for cell in row:
            cell = cell.strip().strip('"')
            if cell and cell.upper() not in {"RUN_ID", "RUN_IDS"}:
                values.append(cell)
    if not values:
        raise ValueError(f"No run IDs found in {path}")
    return values


def collect_metrics_files(input_directory: Path, source: str, run_ids: list[str], metrics_file: Path | None = None) -> list[Path]:
    if metrics_file is not None:
        run_id = run_id_from_filename(metrics_file)
        if run_id not in set(run_ids):
            raise ValueError(f"Input metrics file run ID is not present in RUN_IDs.csv: {run_id}")
        return [metrics_file]

    source_directory = input_directory / source
    if not source_directory.is_dir():
        raise FileNotFoundError(f"Input directory not found: {source_directory}")

    all_files = {run_id_from_filename(path): path for path in sorted(source_directory.glob("*.tsv"))}
    missing = [run_id for run_id in run_ids if run_id not in all_files]
    if missing:
        raise FileNotFoundError(f"Missing {source} metrics files for RUN_IDs: {', '.join(missing)}")
    return [all_files[run_id] for run_id in run_ids]


def cast_all_columns_to_string(df: pl.DataFrame) -> pl.DataFrame:
    return df.select([pl.col(column).cast(pl.Utf8).alias(column) for column in df.columns])


def concat_diagonal(frames: list[pl.DataFrame]) -> pl.DataFrame:
    if not frames:
        return pl.DataFrame()
    frames = [cast_all_columns_to_string(frame) for frame in frames]
    return pl.concat(frames, how="diagonal").fill_null("")


def ordered_columns(frames: list[pl.DataFrame]) -> list[str]:
    columns = []
    for col in LOCALAPP_CANONICAL_COLUMNS:
        if any(col in frame.columns for frame in frames):
            columns.append(col)
    for frame in frames:
        for col in frame.columns:
            if col not in columns:
                columns.append(col)
    return columns


def add_missing_columns(df: pl.DataFrame, columns: list[str], fill_value: str) -> pl.DataFrame:
    expressions = []
    existing = set(df.columns)
    for col in columns:
        if col not in existing:
            expressions.append(pl.lit(fill_value, dtype=pl.Utf8).alias(col))
    if not expressions:
        return df
    return df.with_columns(expressions)


def rename_existing_columns(df: pl.DataFrame, columns: list[str], suffix: str) -> pl.DataFrame:
    rename_map = {col: f"{col}_{suffix}" for col in columns if col in df.columns}
    return df.rename(rename_map)


def build_combined_master(input_df, common_columns, dragen_unique_columns, localapp_unique_columns, source):
    df = input_df.clone()
    df.columns = make_unique_columns(df.columns)
    df = add_missing_columns(df, common_columns, "")
    df = add_missing_columns(df, dragen_unique_columns, "NA" if source == "localapp" else "")
    df = add_missing_columns(df, localapp_unique_columns, "NA" if source == "dragen" else "")
    df = rename_existing_columns(df, dragen_unique_columns, "dragen")
    df = rename_existing_columns(df, localapp_unique_columns, "localapp")
    dragen_output_columns = [f"{col}_dragen" for col in dragen_unique_columns]
    localapp_output_columns = [f"{col}_localapp" for col in localapp_unique_columns]
    required_output_columns = common_columns + dragen_output_columns + localapp_output_columns
    df = add_missing_columns(df, required_output_columns, "NA")
    return df.select(required_output_columns)


def first_value(df: pl.DataFrame, column: str) -> str:
    if column not in df.columns or df.height == 0:
        return "NA"
    for value in df[column].to_list():
        value = clean_cell(value)
        if value != "":
            return value
    return "NA"


def write_per_run_metrics(parsed_frames: list[pl.DataFrame], run_ids: list[str], output_directory: Path) -> list[Path]:
    intermediate_dir = output_directory / "intermediate_metrics_files"
    intermediate_dir.mkdir(parents=True, exist_ok=True)
    paths = []
    for parsed, run_id in zip(parsed_frames, run_ids):
        path = intermediate_dir / f"{run_id}_metrics.tsv"
        parsed.write_csv(path, separator="\t")
        paths.append(path)
    return paths


def write_master(master: pl.DataFrame, output_directory: Path) -> Path:
    intermediate_dir = output_directory / "intermediate_metrics_files"
    intermediate_dir.mkdir(parents=True, exist_ok=True)
    path = intermediate_dir / "master_metrics_table.tsv"
    master.write_csv(path, separator="\t")
    return path


def write_joint_qc_file(parsed_frames: list[pl.DataFrame], run_ids: list[str], output_directory: Path) -> Path:
    intermediate_dir = output_directory / "intermediate_metrics_files"
    intermediate_dir.mkdir(parents=True, exist_ok=True)
    records = []
    for index, (parsed, run_id) in enumerate(zip(parsed_frames, run_ids), start=1):
        records.append(
            {
                "RUN_ID": run_id,
                "PCT_PF_READS": first_value(parsed, "PCT_PF_READS"),
                "PCT_Q30_R1": first_value(parsed, "PCT_Q30_R1"),
                "PCT_Q30_R2": first_value(parsed, "PCT_Q30_R2"),
                "CLUSTER_DENSITY": "NA",
                "ESTIMATED_YIELD": "NA",
                "CLUSTERS_PASSING_FILTER": "NA",
                "RUN_NUMBER": str(index),
            }
        )
    path = intermediate_dir / "joint_sequencing_QC_file.tsv"
    pl.DataFrame(records).write_csv(path, separator="\t")
    return path


def run_plotting_script(script_path: Path, output_directory: Path, master_path: Path, run_metrics_paths: list[Path], joint_qc_path: Path, run_ids: list[str], create_plots: str):
    if create_plots == "False":
        return
    if not script_path.is_file():
        raise FileNotFoundError(f"R plotting script not found: {script_path}")
    pdf_path = output_directory / "TSO500_run_metrics.pdf"
    highlighted_run = run_ids[-1]
    command = [
        "Rscript",
        str(script_path),
        str(len(run_ids)),
        str(joint_qc_path),
        str(pdf_path),
        str(master_path),
        highlighted_run,
        create_plots,
    ] + [str(path) for path in run_metrics_paths] + run_ids
    subprocess.run(command, check=True)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", choices=["localapp", "dragen"], required=True)
    parser.add_argument("--metrics-file", default=None, type=Path)
    parser.add_argument("--input-directory", default=Path("in"), type=Path)
    parser.add_argument("--output-directory", default=Path("out"), type=Path)
    parser.add_argument("--run-id-file", default=None, type=Path)
    parser.add_argument("--r-script", default=None, type=Path)
    parser.add_argument("--create-plots", default="True", choices=["True", "False"])
    parser.add_argument("--no-duplicate-disambiguation", action="store_true")
    return parser.parse_args()


def parse_reference_frames(input_directory: Path, run_ids: list[str], disambiguate: bool):
    localapp_files = collect_metrics_files(input_directory, "localapp", run_ids)
    dragen_files = collect_metrics_files(input_directory, "dragen", run_ids)
    localapp_frames = [parse_metrics_output(path, run_id, disambiguate) for path, run_id in zip(localapp_files, run_ids)]
    dragen_frames = [parse_metrics_output(path, run_id, disambiguate) for path, run_id in zip(dragen_files, run_ids)]
    return localapp_frames, dragen_frames


def main():
    args = parse_args()
    script_directory = Path(__file__).resolve().parent
    r_script = args.r_script or script_directory / "plot_run_metrics.R"
    run_id_file = args.run_id_file or args.input_directory / "RUN_IDs.csv"
    run_ids = read_run_ids(run_id_file)
    disambiguate = not args.no_duplicate_disambiguation

    localapp_reference_frames, dragen_reference_frames = parse_reference_frames(args.input_directory, run_ids, disambiguate)
    localapp_columns = ordered_columns(localapp_reference_frames)
    dragen_columns = ordered_columns(dragen_reference_frames)
    common_columns = [col for col in localapp_columns if col in dragen_columns]
    dragen_unique_columns = [col for col in dragen_columns if col not in common_columns]
    localapp_unique_columns = [col for col in localapp_columns if col not in common_columns]

    selected_files = collect_metrics_files(args.input_directory, args.source, run_ids, args.metrics_file)
    selected_run_ids = [run_id_from_filename(path) for path in selected_files]
    if args.metrics_file is None:
        selected_run_ids = run_ids

    selected_frames = [parse_metrics_output(path, run_id, disambiguate) for path, run_id in zip(selected_files, selected_run_ids)]
    selected_combined = concat_diagonal(selected_frames)
    master = build_combined_master(selected_combined, common_columns, dragen_unique_columns, localapp_unique_columns, args.source)

    run_metrics_paths = write_per_run_metrics(selected_frames, selected_run_ids, args.output_directory)
    master_path = write_master(master, args.output_directory)
    joint_qc_path = write_joint_qc_file(selected_frames, selected_run_ids, args.output_directory)
    # run_plotting_script(r_script, args.output_directory, master_path, run_metrics_paths, joint_qc_path, selected_run_ids, args.create_plots)

    print(f"source: {args.source}")
    print(f"runs: {len(selected_run_ids)}")
    print(f"common columns: {len(common_columns)}")
    print(f"dragen-only columns: {len(dragen_unique_columns)}")
    print(f"localapp-only columns: {len(localapp_unique_columns)}")
    print(f"wrote: {master_path}")
    print(f"wrote: {joint_qc_path}")
    if args.create_plots == "True":
        print(f"wrote: {args.output_directory / 'TSO500_run_metrics.pdf'}")


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise
