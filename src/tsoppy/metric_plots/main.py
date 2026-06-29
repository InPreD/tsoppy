from __future__ import annotations

import re
from collections import OrderedDict
from pathlib import Path

import polars as pl

from tsoppy.general.file_parser import Parse_section_tsv


METRIC_COL = "Metric (UOM)"
LSL_COL = "LSL Guideline"
USL_COL = "USL Guideline"
VALUE_COL = "Value"

NON_SAMPLE_COLUMNS = {METRIC_COL, LSL_COL, USL_COL, VALUE_COL, "-", ""}
MISSING_VALUES = {"", "-", "NA", "N/A", "nan", "None", None}

METADATA_COLUMNS = [
    "SAMPLE_ID",
    "RUN",
    "WORKFLOW_TYPE",
    "WORKFLOW_VERSION",
    "RECORD_TYPE",
]


class MetricPlots:
    def __init__(
        self,
        input_directory: Path,
        run_id_file: Path,
        output_directory: Path,
        create_plots: bool = True,
    ):
        self.input_directory = Path(input_directory)
        self.run_id_file = Path(run_id_file)
        self.output_directory = Path(output_directory)
        self.create_plots = create_plots
        self.intermediate_directory = (
            self.output_directory / "intermediate_metrics_files"
        )

    def run(self):
        run_ids = self._read_run_ids()
        metrics_files = self._find_metrics_files(run_ids)
        parsed_frames = [self._parse_metrics_file(
            path) for path in metrics_files]
        master = self._concat_frames(parsed_frames)

        self.intermediate_directory.mkdir(parents=True, exist_ok=True)

        master_path = self.intermediate_directory / "master_metrics_table.tsv"
        joint_qc_path = self.intermediate_directory / "joint_sequencing_QC_file.tsv"

        self._write_master(master, master_path)
        self._write_joint_sequencing_qc_file(
            parsed_frames, metrics_files, joint_qc_path
        )

        if self.create_plots:
            self._run_plotting_script(
                master_path, joint_qc_path, metrics_files)

    def _read_run_ids(self) -> list[str]:
        if not self.run_id_file.is_file():
            raise FileNotFoundError(
                f"RUN ID file not found: {self.run_id_file}")

        run_ids = []
        with self.run_id_file.open() as handle:
            for line in handle:
                run_id = line.strip()
                if run_id and not run_id.startswith("#"):
                    run_ids.append(run_id)

        if not run_ids:
            raise ValueError(f"No RUN IDs found in {self.run_id_file}")

        return run_ids

    def _find_metrics_files(self, run_ids: list[str]) -> list[Path]:
        files = []

        for run_id in run_ids:
            matches = []
            for subdir in ["dragen", "localapp"]:
                directory = self.input_directory / subdir
                print(directory)
                if directory.is_dir():
                    matches.extend(
                        sorted(directory.glob(f"{run_id}*MetricsOutput*.tsv"))
                    )

            if not matches:
                raise FileNotFoundError(
                    f"No MetricsOutput.tsv found for RUN ID: {run_id}"
                )

            files.extend(matches)

        return files

    def _parse_metrics_file(self, path: Path) -> pl.DataFrame:
        headers, sections = Parse_section_tsv(str(path), ["Header"])
        workflow_type, workflow_version = self._detect_workflow(
            headers, sections, path)
        run_id = self._run_id_from_filename(path)

        sample_records: OrderedDict[str, OrderedDict[str, str]] = OrderedDict()
        run_metrics: OrderedDict[str, str] = OrderedDict()
        lsl_metrics: OrderedDict[str, str] = OrderedDict()
        usl_metrics: OrderedDict[str, str] = OrderedDict()

        for section_name, df in sections.items():
            if df.is_empty() or section_name == "Header":
                continue

            metric_col = METRIC_COL if METRIC_COL in df.columns else df.columns[0]
            samples = self._sample_columns(df)

            if VALUE_COL in df.columns and not samples:
                for row in df.iter_rows(named=True):
                    metric = self._metric_name_for_section(
                        row.get(metric_col, ""), section_name
                    )
                    if metric:
                        self._merge_value(run_metrics, metric,
                                          row.get(VALUE_COL, ""))
                continue

            for row in df.iter_rows(named=True):
                metric = self._metric_name_for_section(
                    row.get(metric_col, ""), section_name
                )
                if not metric:
                    continue

                if LSL_COL in df.columns:
                    self._merge_value(lsl_metrics, metric,
                                      row.get(LSL_COL, "NA"))

                if USL_COL in df.columns:
                    self._merge_value(usl_metrics, metric,
                                      row.get(USL_COL, "NA"))

                for sample_id in samples:
                    sample_id = self._clean(sample_id)
                    if not sample_id or sample_id in NON_SAMPLE_COLUMNS:
                        continue

                    sample_records.setdefault(sample_id, OrderedDict())
                    self._merge_value(
                        sample_records[sample_id], metric, row.get(
                            sample_id, "")
                    )

        records = []

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
                record[metric] = self._clean(value) or "NA"

            records.append(record)

        for sample_id, metrics in sample_records.items():
            record = OrderedDict()
            record["SAMPLE_ID"] = sample_id
            record["RUN"] = run_id
            record["WORKFLOW_TYPE"] = workflow_type
            record["WORKFLOW_VERSION"] = workflow_version
            record["RECORD_TYPE"] = self._infer_record_type(sample_id, metrics)

            for metric, value in run_metrics.items():
                record[metric] = value

            for metric, value in metrics.items():
                record[metric] = value

            records.append(record)

        return pl.DataFrame(records, infer_schema_length=None)

    def _write_master(self, master: pl.DataFrame, output_path: Path) -> None:
        master.select(self._order_columns(master)).write_csv(
            output_path, separator="\t"
        )

    def _write_joint_sequencing_qc_file(
        self,
        parsed_frames: list[pl.DataFrame],
        metrics_files: list[Path],
        output_path: Path,
    ) -> None:
        records = []

        for run_number, (frame, path) in enumerate(
            zip(parsed_frames, metrics_files), start=1
        ):
            records.append(
                {
                    "RUN_ID": self._run_id_from_filename(path),
                    "PCT_PF_READS": self._first_sample_value(frame, "PCT_PF_READS"),
                    "PCT_Q30_R1": self._first_sample_value(frame, "PCT_Q30_R1"),
                    "PCT_Q30_R2": self._first_sample_value(frame, "PCT_Q30_R2"),
                    "CLUSTER_DENSITY": "0",
                    "ESTIMATED_YIELD": "0",
                    "CLUSTERS_PASSING_FILTER": "0",
                    "RUN_NUMBER": str(run_number),
                }
            )

        pl.DataFrame(records).select(
            [
                "RUN_ID",
                "PCT_PF_READS",
                "PCT_Q30_R1",
                "PCT_Q30_R2",
                "CLUSTER_DENSITY",
                "ESTIMATED_YIELD",
                "CLUSTERS_PASSING_FILTER",
                "RUN_NUMBER",
            ]
        ).write_csv(output_path, separator="\t")

    def _run_plotting_script(
        self,
        master_path: Path,
        joint_qc_path: Path,
        metrics_files: list[Path],
    ) -> None:
        r_script = Path("plot_run_metrics.R")
        if not r_script.is_file():
            return

    @staticmethod
    def _clean(value: object) -> str:
        if value is None:
            return ""
        return str(value).strip()

    @staticmethod
    def _is_missing(value: object) -> bool:
        if value is None:
            return True
        return str(value).strip() in MISSING_VALUES

    @staticmethod
    def _normalize_metric_name(metric: object) -> str:
        metric = "" if metric is None else str(metric).strip()
        metric = re.sub(r"\s*\([^)]*\)\s*$", "", metric)
        metric = metric.replace("%", "PCT")
        metric = re.sub(r"[^0-9A-Za-z]+", "_", metric)
        metric = re.sub(r"_+", "_", metric).strip("_")
        return metric

    def _metric_name_for_section(self, metric: object, section: str) -> str:
        metric = self._normalize_metric_name(metric)
        section_upper = self._normalize_metric_name(section).upper()

        if "DNA" in section_upper and metric and not metric.startswith("DNA_"):
            return f"DNA_{metric}"

        if "RNA" in section_upper and metric and not metric.startswith("RNA_"):
            return f"RNA_{metric}"

        return metric

    def _sample_columns(self, df: pl.DataFrame) -> list[str]:
        return [col for col in df.columns if self._clean(col) not in NON_SAMPLE_COLUMNS]

    def _merge_value(
        self, record: OrderedDict[str, str], key: str, value: object
    ) -> None:
        value = self._clean(value)

        if key not in record:
            record[key] = value
            return

        if self._is_missing(record[key]) and not self._is_missing(value):
            record[key] = value

    def _infer_record_type(self, sample_id: str, record: OrderedDict[str, str]) -> str:
        sample = sample_id.upper()

        if (
            "_D" in sample
            or "-D" in sample
            or sample.startswith("TVD")
            or sample.startswith("DNA")
        ):
            return "DNA_SAMPLE"

        if (
            "_R" in sample
            or "-R" in sample
            or sample.startswith("TVR")
            or sample.startswith("RNA")
        ):
            return "RNA_SAMPLE"

        dna_count = sum(
            1
            for key, value in record.items()
            if key.startswith("DNA_") and not self._is_missing(value)
        )
        rna_count = sum(
            1
            for key, value in record.items()
            if key.startswith("RNA_") and not self._is_missing(value)
        )

        if dna_count > rna_count:
            return "DNA_SAMPLE"

        if rna_count > dna_count:
            return "RNA_SAMPLE"

        return "SAMPLE"

    def _detect_workflow(
        self,
        headers: list[str],
        sections: dict[str, pl.DataFrame],
        path: Path,
    ) -> tuple[str, str]:
        header_text = " ".join(headers).lower()
        workflow_type = (
            "dragen"
            if "dragen" in header_text or "dragen" in path.name.lower()
            else "localapp"
        )
        workflow_version = "Unknown"

        header_df = sections.get("Header")
        if header_df is not None and not header_df.is_empty():
            if "Workflow Version" in header_df.columns:
                values = [
                    self._clean(value)
                    for value in header_df["Workflow Version"].to_list()
                    if self._clean(value)
                ]
                if values:
                    workflow_version = values[0]

        if workflow_version != "Unknown" and "dragen" in workflow_version.lower():
            workflow_type = "dragen"

        return workflow_type, workflow_version

    @staticmethod
    def _run_id_from_filename(path: Path) -> str:
        name = path.name
        name = re.sub(
            r"_MetricsOutput_(Localapp|LocalApp|Dragen|DRAGEN)\.tsv$", "", name
        )
        name = re.sub(r"_MetricsOutput\.tsv$", "", name)
        return name

    def _concat_frames(self, frames: list[pl.DataFrame]) -> pl.DataFrame:
        if not frames:
            return pl.DataFrame()

        frames = [self._cast_all_string(frame) for frame in frames]
        return pl.concat(frames, how="diagonal").fill_null("NA")

    @staticmethod
    def _cast_all_string(df: pl.DataFrame) -> pl.DataFrame:
        string_type = pl.String if hasattr(pl, "String") else pl.Utf8
        return df.select(
            [pl.col(col).cast(string_type).alias(col) for col in df.columns]
        )

    @staticmethod
    def _order_columns(df: pl.DataFrame) -> list[str]:
        metadata = [col for col in METADATA_COLUMNS if col in df.columns]

        run_metrics = sorted(
            col
            for col in df.columns
            if col not in metadata
            and not col.startswith("DNA_")
            and not col.startswith("RNA_")
        )

        dna_metrics = sorted(
            col for col in df.columns if col.startswith("DNA_"))
        rna_metrics = sorted(
            col for col in df.columns if col.startswith("RNA_"))

        return metadata + run_metrics + dna_metrics + rna_metrics

    def _first_sample_value(self, frame: pl.DataFrame, column: str) -> str:
        if column not in frame.columns:
            return "0"

        sample_rows = frame.filter(
            ~pl.col("RECORD_TYPE").is_in(
                ["LOWER_THRESHOLD", "UPPER_THRESHOLD"])
        )

        for value in sample_rows[column].to_list():
            value = self._clean(value)
            if value and value != "NA":
                return value

        return "0"
