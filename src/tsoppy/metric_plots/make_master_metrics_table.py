#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import OrderedDict
from pathlib import Path
from typing import Iterable

import polars as pl

try:
    from file_parser import Parse_section_tsv
except ImportError as error:
    raise ImportError(
        "Could not import Parse_section_tsv. Put file_parser.py in the same directory "
        "as this script, or run this script from that directory."
    ) from error

"""
RUN_IDS = [
    "240809_A02134_0013_BHCGJYDRX5",
    "240906_A02134_0019_BHHGKGDRX5",
    "240927_A02134_0024_BHHGG7DRX5",
    "241115_A02134_0038_BHKT5WDRX5",
    "241004_A02134_0026_BHJ2CFDRX5",
    "250321_A02134_0073_BHTMFHDRX5",
    "250404_A02134_0079_BHTLC7DRX5",
    "250425_A02134_0083_BHYHMWDRX5",
    "250516_A02134_0087_BHLMFHDRX5",
    "250530_A02134_0089_BHYGC3DRX5",
    "250613_A02134_0091_BHYFHJDRX5",
    "250627_A02134_0093_BHYHNYDRX5",
    "250822_A02134_0101_BHYFHLDRX5",
]

"""

RUN_IDS = [
    "240809_A02134_0013_BHCGJYDRX5",
    "250822_A02134_0101_BHYFHLDRX5",
]

METRIC_COL = "Metric (UOM)"
LSL_COL = "LSL Guideline"
USL_COL = "USL Guideline"
VALUE_COL = "Value"
NON_SAMPLE_COLUMNS = {METRIC_COL, LSL_COL, USL_COL, VALUE_COL, "-", ""}
MISSING_VALUES = {"", "-", "NA", "N/A", "nan", "None", None}
METADATA_COLUMNS = ["SAMPLE_ID", "RUN", "WORKFLOW_TYPE", "WORKFLOW_VERSION", "RECORD_TYPE"]

def read_run_ids(run_id_file: Path) -> list[str]:
    if not run_id_file.is_file():
        raise FileNotFoundError(
            f"RUN ID file does not exist: {run_id_file}"
        )

    run_ids = []

    with open(run_id_file, "r") as handle:
        for line in handle:
            run_id = line.strip()

            if not run_id:
                continue

            if run_id.startswith("#"):
                continue

            run_ids.append(run_id)

    if not run_ids:
        raise ValueError(
            f"No RUN IDs found in {run_id_file}"
        )

    return run_ids

def clean(value: object) -> str:
    if value is None:
        return ""
    return str(value).strip()


def is_missing(value: object) -> bool:
    return clean(value) in MISSING_VALUES


def normalize_metric_name(metric: object) -> str:
    metric = clean(metric)
    metric = re.sub(r"\s*\([^)]*\)\s*$", "", metric)
    metric = metric.replace("%", "PCT")
    metric = re.sub(r"[^0-9A-Za-z]+", "_", metric)
    metric = re.sub(r"_+", "_", metric).strip("_")
    return metric


def normalize_section_name(section: str) -> str:
    return normalize_metric_name(section).upper()


def metric_prefix(section: str) -> str:
    section = normalize_section_name(section)
    if "DNA" in section:
        return "DNA"
    if "RNA" in section:
        return "RNA"
    return ""


def metric_name_for_section(metric: object, section: str) -> str:
    metric = normalize_metric_name(metric)
    prefix = metric_prefix(section)
    if prefix and metric and not metric.startswith(f"{prefix}_"):
        return f"{prefix}_{metric}"
    return metric


def run_id_from_filename(path: Path) -> str:
    name = path.name
    name = re.sub(r"_MetricsOutput_(Localapp|LocalApp|Dragen|DRAGEN)\.tsv$", "", name)
    name = re.sub(r"_MetricsOutput\.tsv$", "", name)
    return name


def make_unique_headers(headers: Iterable[object]) -> list[str]:
    seen: dict[str, int] = {}
    out: list[str] = []
    for header in headers:
        header = clean(header) or "-"
        if header not in seen:
            seen[header] = 0
            out.append(header)
        else:
            seen[header] += 1
            out.append(f"{header}_{seen[header]}")
    return out


def parse_section_tsv_safe(path: str, key_value_sections: list[str]) -> tuple[list[str], dict[str, pl.DataFrame]]:
    df = pl.read_csv(path, separator="\t", has_header=False)
    section_idx = get_section_idx_safe(df)
    headers: list[str] = []

    if section_idx and section_idx[0][1] != 1:
        for row in df.head(section_idx[0][1] - 1).rows():
            for value in row:
                if value is not None:
                    headers.append(clean(value))

    sections: dict[str, pl.DataFrame] = {}
    for section_name, start, length in section_idx:
        df_slice = df.slice(start, length)

        df_slice = df_slice.select(
            [pl.col(col) for col in df_slice.columns if df_slice[col].null_count() != df_slice.height]
        )

        string_type = pl.String if hasattr(pl, "String") else pl.Utf8
        df_slice = df_slice.with_columns(pl.all().cast(string_type).fill_null("-"))

        if section_name in key_value_sections and df_slice.width == 2:
            df_slice = df_slice.transpose()

        raw_header = list(df_slice.row(0))
        unique_header = make_unique_headers(raw_header)
        rename_map = {old: new for old, new in zip(df_slice.columns, unique_header)}
        sections[section_name] = df_slice.rename(rename_map).slice(1)

    return headers, sections


def get_section_idx_safe(df: pl.DataFrame) -> list[tuple[str, int, int]]:
    section = ""
    section_start = 0
    section_idx: list[tuple[str, int, int]] = []
    rows = df.with_row_index().iter_rows() if hasattr(df, "with_row_index") else df.with_row_count().iter_rows()

    for row in rows:
        row_index = row[0]
        first_col = row[1] if len(row) > 1 else None

        if first_col:
            match = re.search(r"^\[(?P<section>.*)\]$", str(first_col))
            if match:
                section_start = row_index + 1
                section = match.group("section")

        if all(item is None for item in row[1:]):
            section_length = row_index - section_start
            if section_start > 0 and section_length > 0:
                section_idx.append((section, section_start, section_length))
                section_start = 0

        if row_index == len(df) - 1:
            section_length = row_index - section_start + 1
            if section_start > 0 and section_length > 0:
                section_idx.append((section, section_start, section_length))

    return section_idx


def parse_sections(path: Path) -> tuple[list[str], dict[str, pl.DataFrame]]:
    return Parse_section_tsv(str(path), ["Header"])
"""
    try:
        return Parse_section_tsv(str(path), ["Header"])
    except Exception:
        return parse_section_tsv_safe(str(path), ["Header"])
"""

def detect_workflow(headers: list[str], sections: dict[str, pl.DataFrame], path: Path) -> tuple[str, str]:
    header_text = " ".join(headers).lower()
    workflow_type = "dragen" if "dragen" in header_text or "dragen" in path.name.lower() else "localapp"

    workflow_version = "Unknown"
    header_df = sections.get("Header")
    if header_df is not None and not header_df.is_empty():
        if "Workflow Version" in header_df.columns:
            values = [clean(v) for v in header_df["Workflow Version"].to_list() if clean(v)]
            if values:
                workflow_version = values[0]

        if workflow_version == "Unknown":
            for row in header_df.iter_rows(named=True):
                values = list(row.values())
                for index, value in enumerate(values):
                    value = clean(value)
                    if value.lower() == "workflow version" and index + 1 < len(values):
                        workflow_version = clean(values[index + 1])
                    elif "workflow version" in value.lower():
                        parts = re.split(r"[:=]", value, maxsplit=1)
                        if len(parts) == 2:
                            workflow_version = parts[1].strip()

    if workflow_version != "Unknown" and "dragen" in workflow_version.lower():
        workflow_type = "dragen"

    return workflow_type, workflow_version


def sample_columns(df: pl.DataFrame) -> list[str]:
    return [col for col in df.columns if clean(col) not in NON_SAMPLE_COLUMNS]


def merge_value(record: OrderedDict[str, str], key: str, value: object) -> None:
    value = clean(value)
    if key not in record:
        record[key] = value
        return
    if is_missing(record[key]) and not is_missing(value):
        record[key] = value


def infer_record_type(sample_id: str, record: OrderedDict[str, str]) -> str:
    sample = sample_id.upper()
    if "_D" in sample or "-D" in sample or sample.startswith("DNA"):
        return "DNA_SAMPLE"
    if "_R" in sample or "-R" in sample or sample.startswith("RNA"):
        return "RNA_SAMPLE"

    dna_count = sum(1 for k, v in record.items() if k.startswith("DNA_") and not is_missing(v))
    rna_count = sum(1 for k, v in record.items() if k.startswith("RNA_") and not is_missing(v))
    if dna_count > rna_count:
        return "DNA_SAMPLE"
    if rna_count > dna_count:
        return "RNA_SAMPLE"
    return "SAMPLE"


def parse_metrics_file(path: Path) -> pl.DataFrame:
    headers, sections = parse_sections(path)
    workflow_type, workflow_version = detect_workflow(headers, sections, path)
    run_id = run_id_from_filename(path)

    sample_records: OrderedDict[str, OrderedDict[str, str]] = OrderedDict()
    run_metrics: OrderedDict[str, str] = OrderedDict()
    lsl_metrics: OrderedDict[str, str] = OrderedDict()
    usl_metrics: OrderedDict[str, str] = OrderedDict()

    for section_name, df in sections.items():
        if df.is_empty() or section_name == "Header":
            continue
        
        metric_col = METRIC_COL if METRIC_COL in df.columns else df.columns[0]
        samples = sample_columns(df)

        if VALUE_COL in df.columns and not samples:
            for row in df.iter_rows(named=True):
                metric = metric_name_for_section(row.get(metric_col, ""), section_name)
                if metric:
                    merge_value(run_metrics, metric, row.get(VALUE_COL, ""))
            continue

        for row in df.iter_rows(named=True):
            metric = metric_name_for_section(row.get(metric_col, ""), section_name)
            if not metric:
                continue

            if LSL_COL in df.columns:
                merge_value(lsl_metrics, metric, row.get(LSL_COL, "NA"))
            if USL_COL in df.columns:
                merge_value(usl_metrics, metric, row.get(USL_COL, "NA"))

            for sample_id in samples:
                sample_id = clean(sample_id)
                if not sample_id or sample_id in NON_SAMPLE_COLUMNS:
                    continue
                sample_records.setdefault(sample_id, OrderedDict())
                merge_value(sample_records[sample_id], metric, row.get(sample_id, ""))

    records: list[OrderedDict[str, str]] = []

    for sample_id, record_type, metric_values in [
        ("LSL_Guideline", "LOWER_THRESHOLD", lsl_metrics),
        ("USL_Guideline", "UPPER_THRESHOLD", usl_metrics),
    ]:
        record = OrderedDict()
        record["SAMPLE_ID"] = sample_id
        record["RUN"] = run_id
        record["WORKFLOW_TYPE"] = workflow_type
        record["WORKFLOW_VERSION"] = workflow_version
        record["RECORD_TYPE"] = record_type
        for metric in run_metrics:
            record[metric] = "NA"
        for metric, value in metric_values.items():
            record[metric] = clean(value) or "NA"
        records.append(record)

    for sample_id, metrics in sample_records.items():
        record = OrderedDict()
        record["SAMPLE_ID"] = sample_id
        record["RUN"] = run_id
        record["WORKFLOW_TYPE"] = workflow_type
        record["WORKFLOW_VERSION"] = workflow_version
        record["RECORD_TYPE"] = infer_record_type(sample_id, metrics)
        for metric, value in run_metrics.items():
            record[metric] = value
        for metric, value in metrics.items():
            record[metric] = value
        records.append(record)

    return pl.DataFrame(records, infer_schema_length=None)


def find_metrics_files(input_dir: Path, run_ids: list[str]) -> list[Path]:
    candidates: list[Path] = []
    for run_id in run_ids:
        matches = sorted((input_dir / "dragen").glob(f"{run_id}*MetricsOutput*.tsv"))
        matches += sorted((input_dir / "localapp").glob(f"{run_id}*MetricsOutput*.tsv"))
        if not matches:
            raise FileNotFoundError(f"No MetricsOutput.tsv found for RUN_ID: {run_id}")
        candidates.extend(matches)
    return candidates


def cast_all_string(df: pl.DataFrame) -> pl.DataFrame:
    string_type = pl.String if hasattr(pl, "String") else pl.Utf8
    return df.select([pl.col(col).cast(string_type).alias(col) for col in df.columns])


def concat_frames(frames: list[pl.DataFrame]) -> pl.DataFrame:
    frames = [cast_all_string(frame) for frame in frames]
    return pl.concat(frames, how="diagonal").fill_null("NA")


def order_columns(df: pl.DataFrame) -> list[str]:
    metadata = [col for col in METADATA_COLUMNS if col in df.columns]
    run_metrics = sorted(
        col for col in df.columns
        if col not in metadata and not col.startswith("DNA_") and not col.startswith("RNA_")
    )
    dna_metrics = sorted(col for col in df.columns if col.startswith("DNA_"))
    rna_metrics = sorted(col for col in df.columns if col.startswith("RNA_"))
    return metadata + run_metrics + dna_metrics + rna_metrics


def write_joint_metrics_header(path: Path) -> None:
    path.write_text(
        "# Joint metrics table generated from MetricsOutput.tsv files parsed with file_parser.py / Parse_section_tsv.\n"
        "# WORKFLOW_TYPE and WORKFLOW_VERSION identify whether each row came from DRAGEN or LocalApp.\n"
        "# RECORD_TYPE identifies DNA_SAMPLE, RNA_SAMPLE, LOWER_THRESHOLD, or UPPER_THRESHOLD rows.\n"
        "# DNA-specific metric columns use the DNA_ prefix. RNA-specific metric columns use the RNA_ prefix. Metadata columns do not use DNA_ or RNA_ prefixes.\n"
        "# Interpret metrics cautiously when multiple workflows contribute data because similarly named metrics can be generated or formatted differently across workflows.\n"
    )


def write_master(master: pl.DataFrame, output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    ordered = master.select(order_columns(master))
    ordered.write_csv(output_path, separator="\t")


def write_joint_metrics_file(master: pl.DataFrame, output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    write_joint_metrics_header(output_path)
    with output_path.open("a") as handle:
        handle.write(master.select(order_columns(master)).write_csv(separator="\t"))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-directory", type=Path, default=Path("in"))
    parser.add_argument("--output", type=Path, default=Path("out/intermediate_metrics_files/master_metrics_table.tsv"))
    parser.add_argument("--joint-output", type=Path, default=Path("out/intermediate_metrics_files/joint_metrics_file.tsv"))
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    files = find_metrics_files(args.input_directory, RUN_IDS)
    frames = [parse_metrics_file(path) for path in files]
    master = concat_frames(frames)

    write_master(master, args.output)
    write_joint_metrics_file(master, args.joint_output)

    print(f"files parsed: {len(files)}")
    print(f"rows: {master.height}")
    print(f"columns: {master.width}")
    print(f"wrote: {args.output}")
    print(f"wrote: {args.joint_output}")


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise
