"""Create standardized TSO500 metrics tables."""

from __future__ import annotations

import logging
from pathlib import Path
from glob import glob
import polars as pl

from tsoppy.general.classes import MetricsOutputTsv, WorkflowOutput


logger = logging.getLogger(__name__)


class MetricPlots:
    """Create metrics tables and prepare data for QC plotting.

    Attributes:
        config_yaml: Workflow output configuration path.
        inpred_nomenclature: InPreD nomenclature configuration path.
        input_glob: Glob pattern matching workflow output roots.
        input_roots: Workflow roots resolved from input_glob.
        run_ids: Run IDs selected for processing.
    """

    METRIC_COL = "Metric (UOM)"
    LSL_COL = "LSL Guideline"
    USL_COL = "USL Guideline"
    VALUE_COL = "Value"

    LSL_SAMPLE_ID = "LSL_Guideline"
    USL_SAMPLE_ID = "USL_Guideline"
    RUN_VALUE_ID = "__RUN_VALUE__"

    LOWER_THRESHOLD = "LOWER_THRESHOLD"
    UPPER_THRESHOLD = "UPPER_THRESHOLD"
    DNA_SAMPLE = "DNA_SAMPLE"
    RNA_SAMPLE = "RNA_SAMPLE"
    UNKNOWN_SAMPLE = "SAMPLE"

    RUN_INDEX = "RUN_INDEX"

    IGNORED_SECTIONS = {
        "Header",
        "Notes",
    }

    MISSING_VALUES = [
        "",
        "-",
        "NA",
        "N/A",
        "nan",
        "None",
    ]

    METADATA_COLUMNS = [
        RUN_INDEX,
        "SAMPLE_ID",
        "RUN",
        "WORKFLOW_TYPE",
        "WORKFLOW_VERSION",
        "RECORD_TYPE",
    ]

    JOINT_QC_METRICS = [
        "PCT_PF_READS",
        "PCT_Q30_R1",
        "PCT_Q30_R2",
        "CLUSTER_DENSITY",
        "ESTIMATED_YIELD",
        "CLUSTERS_PASSING_FILTER",
    ]

    JOINT_QC_COLUMNS = [
        "RUN_INDEX",
        "RUN_ID",
        "WORKFLOW_TYPE",
        "WORKFLOW_VERSION",
        "PCT_PF_READS",
        "PCT_Q30_R1",
        "PCT_Q30_R2",
        "CLUSTER_DENSITY",
        "ESTIMATED_YIELD",
        "CLUSTERS_PASSING_FILTER",
    ]

    def __init__(
        self,
        config_yaml: Path,
        inpred_nomenclature: Path,
        input_glob: str,
        run_ids: str | list[str] | None = None,
        run_id_file: Path | None = None,
    ):
        """Initialize metric plot processing.

        Args:
            config_yaml: Workflow configuration file.
            inpred_nomenclature: InPreD nomenclature file.
            input_glob: Glob matching workflow output roots.
            workdir: Output directory.
            run_ids: Comma-separated string or list of run IDs.
            run_id_file: File containing run IDs.
        """

        # although typer checks for mutual exclusivity, we also check here to ensure that the class is used correctly from other contexts
        if run_ids is not None and run_id_file is not None:
            raise ValueError("run_ids and run_id_file are mutually exclusive.")

        self.config_yaml = Path(config_yaml)
        self.inpred_nomenclature = Path(inpred_nomenclature)
        self.input_glob = input_glob

        self.run_ids = self._resolve_run_ids(
            run_ids=run_ids,
            run_id_file=run_id_file,
        )

        self.input_roots: list[Path] = []

    def generate_metrics_tables(
        self,
    ) -> tuple[
        pl.DataFrame,
        pl.DataFrame,
    ]:
        """Create and write master and joint QC tables."""
        logger.info(
            "Generating metrics tables for %d run(s).",
            len(self.run_ids),
        )
        metrics_outputs = self._load_metrics_outputs()

        run_frames = [
            self._transform_metrics_output(
                run_id=run_id,
                metrics_output=metrics_output,
            )
            for run_id, metrics_output in metrics_outputs
        ]

        master = self._combine_runs(run_frames)
        master = self._add_run_index(
            master,
            run_column="RUN",
        )
        master = self._finalize_frame(master)

        joint_qc = self._create_joint_qc(master)

        master_path = Path("master_metrics_table.tsv")
        joint_qc_path = Path("joint_sequencing_QC_file.tsv")

        master.write_csv(
            master_path,
            separator="\t",
        )
        joint_qc.write_csv(
            joint_qc_path,
            separator="\t",
        )

        logger.info(
            "Generated %d master rows and %d joint QC rows.",
            master.height,
            joint_qc.height,
        )

        return master, joint_qc

    def select_plot_data(
        self,
        master: pl.DataFrame,
        joint_qc: pl.DataFrame,
        workflow_type: str,
        plot_last_runs: int | None = None,
        plot_run_ids: list[str] | None = None,
    ) -> tuple[pl.DataFrame, pl.DataFrame]:
        """Prepare workflow-specific metrics and joint QC rows for plotting."""

        workflow_frame = master.filter(pl.col("WORKFLOW_TYPE") == workflow_type)

        available_runs = (
            workflow_frame.select("RUN")
            .unique(maintain_order=True)
            .get_column("RUN")
            .to_list()
        )

        if plot_last_runs is not None:
            selected_runs = available_runs[-plot_last_runs:]

        elif plot_run_ids is not None:
            missing_runs = [
                run_id for run_id in plot_run_ids if run_id not in available_runs
            ]

            if missing_runs:
                logger.warning(
                    "The following run IDs are not available for %s "
                    "and will be skipped: %s",
                    workflow_type,
                    ", ".join(missing_runs),
                )

            selected_runs = [
                run_id for run_id in plot_run_ids if run_id in available_runs
            ]

        else:
            selected_runs = available_runs

        logger.info(
            "Selected %d %s run(s) for plotting.",
            len(selected_runs),
            workflow_type,
        )

        plot_frame = workflow_frame.filter(pl.col("RUN").is_in(selected_runs))

        plot_joint_qc = joint_qc.filter(
            (pl.col("WORKFLOW_TYPE") == workflow_type)
            & pl.col("RUN_ID").is_in(selected_runs)
        )

        return plot_frame, plot_joint_qc

    def _resolve_run_ids(
        self,
        run_ids: str | list[str] | None,
        run_id_file: Path | None,
    ) -> list[str]:
        """Read run IDs from CLI values or a run ID file."""

        message = "No run IDs were provided. Use run_ids or run_id_file."

        if run_ids is not None:
            selected_run_ids = self._parse_run_id_input(run_ids)

        elif run_id_file is not None:
            selected_run_ids = self._read_run_id_file(run_id_file)

        else:
            logger.error(message)
            raise ValueError(message)

        unique_run_ids = list(dict.fromkeys(selected_run_ids))

        if not unique_run_ids:
            raise ValueError(message)

        logger.info(
            "Selected %d unique run ID(s).",
            len(unique_run_ids),
        )

        return unique_run_ids

    def _read_run_id_file(
        self,
        run_id_file: Path,
    ) -> list[str]:
        """Read one or more run IDs from a text file."""

        if not run_id_file.is_file():
            message = f"Run ID file does not exist: {run_id_file}."
            logger.error(message)
            raise FileNotFoundError(message)

        run_ids: list[str] = []

        with run_id_file.open(
            encoding="utf-8",
        ) as handle:
            for line in handle:
                cleaned = line.strip()

                if not cleaned or cleaned.startswith("#"):
                    continue

                run_ids.extend(self._parse_run_id_input(cleaned))

        return run_ids

    @staticmethod
    def _parse_run_id_input(
        run_ids: str | list[str],
    ) -> list[str]:
        """Parse comma-separated or list-based run IDs."""
        if isinstance(run_ids, str):
            values = run_ids.split(",")
        else:
            values = [item for value in run_ids for item in value.split(",")]

        return [
            value.strip().strip('"').strip("'")
            for value in values
            if value.strip() and not value.strip().startswith("#")
        ]

    def _load_metrics_outputs(
        self,
    ) -> list[tuple[str, MetricsOutputTsv]]:
        """Load all configured workflow outputs for selected runs."""
        self._validate_inputs()

        outputs: list[tuple[str, MetricsOutputTsv]] = []

        for run_id in self.run_ids:
            run_roots = self._find_run_roots(run_id)

            if not run_roots:
                message = f"No workflow output root found for run {run_id}."
                logger.error(message)
                raise FileNotFoundError(message)

            loaded_for_run = 0

            for root in run_roots:
                try:
                    workflow_output = WorkflowOutput(
                        config_yaml=self.config_yaml,
                        inpred_nomenclature=self.inpred_nomenclature,
                        root_path=root,
                    )

                    metrics_output = MetricsOutputTsv.create(workflow_output)

                except (
                    FileNotFoundError,
                    KeyError,
                    ValueError,
                ) as error:
                    logger.debug(f"Skipping workflow root {root}: {error}")
                    continue

                outputs.append(
                    (
                        run_id,
                        metrics_output,
                    )
                )
                loaded_for_run += 1

                logger.info(
                    f"Loaded {metrics_output.workflow_type} "
                    f"{metrics_output.workflow_version} "
                    f"for run {run_id} from {metrics_output.path}."
                )

            if loaded_for_run == 0:
                message = f"No valid MetricsOutput.tsv found for run {run_id}."
                logger.error(message)
                raise FileNotFoundError(message)

        return outputs

    def _validate_inputs(
        self,
    ) -> None:
        """Validate required inputs and resolve workflow roots."""

        for path, description in [
            (
                self.config_yaml,
                "workflow configuration",
            ),
            (
                self.inpred_nomenclature,
                "InPreD nomenclature",
            ),
        ]:
            if not path.is_file():
                message = f"{description} file does not exist: {path}"
                logger.error(message)
                raise FileNotFoundError(message)

        self.input_roots = list(
            dict.fromkeys(
                Path(match)
                for match in glob(
                    self.input_glob,
                    recursive=True,
                )
                if Path(match).is_dir()
            )
        )

        if not self.input_roots:
            message = (
                f"Input glob did not match any workflow output directories: "
                f"{self.input_glob}"
            )
            logger.error(message)
            raise FileNotFoundError(message)

        logger.info(f"Input glob matched {len(self.input_roots)} workflow root(s).")

    def _find_run_roots(
        self,
        run_id: str,
    ) -> list[Path]:
        """Find glob-matched workflow roots for a selected run ID."""

        return [root for root in self.input_roots if root.name == run_id]

    def _transform_metrics_output(
        self,
        run_id: str,
        metrics_output: MetricsOutputTsv,
    ) -> pl.DataFrame:
        """Transform one parsed MetricsOutput.tsv."""
        sample_frames: list[pl.DataFrame] = []
        threshold_frames: list[pl.DataFrame] = []
        run_metric_frames: list[pl.DataFrame] = []

        for (
            section_name,
            section,
        ) in metrics_output.sections.items():
            if section_name in self.IGNORED_SECTIONS or section.is_empty():
                continue

            standardized = self._standardize_section(
                section_name=section_name,
                section=section,
            )

            sample_columns = self._sample_columns(standardized)

            if sample_columns:
                sample_frames.append(
                    self._section_to_wide(
                        section=standardized,
                        value_columns=(sample_columns),
                    )
                )

            threshold_frames.append(
                self._section_to_wide(
                    section=standardized,
                    value_columns=[
                        self.LSL_COL,
                        self.USL_COL,
                    ],
                    identifier_mapping={
                        self.LSL_COL: (self.LSL_SAMPLE_ID),
                        self.USL_COL: (self.USL_SAMPLE_ID),
                    },
                )
            )

            if self.VALUE_COL in standardized.columns:
                run_metric_frames.append(
                    self._section_to_wide(
                        section=standardized,
                        value_columns=[self.VALUE_COL],
                        identifier_mapping={self.VALUE_COL: (self.RUN_VALUE_ID)},
                    )
                )

        samples = self._merge_wide_frames(sample_frames)
        thresholds = self._merge_wide_frames(threshold_frames)
        run_metrics = self._merge_wide_frames(run_metric_frames)

        if not run_metrics.is_empty():
            if run_metrics.height != 1:
                message = (
                    f"Expected exactly one run-level metrics row for run {run_id}, "
                    f"found {run_metrics.height}."
                )
                logger.error(message)
                raise ValueError(message)

            run_values = run_metrics.drop("SAMPLE_ID")

            if samples.is_empty():
                logger.warning(f"No sample rows found for run {run_id}.")
            else:
                samples = samples.join(
                    run_values,
                    how="cross",
                )

        if not samples.is_empty():
            samples = self._add_record_type(samples)

        if not thresholds.is_empty():
            thresholds = thresholds.with_columns(
                pl.when(pl.col("SAMPLE_ID") == self.LSL_SAMPLE_ID)
                .then(pl.lit(self.LOWER_THRESHOLD))
                .otherwise(pl.lit(self.UPPER_THRESHOLD))
                .alias("RECORD_TYPE")
            )

        output_frames = [
            frame
            for frame in [
                thresholds,
                samples,
            ]
            if not frame.is_empty()
        ]

        if not output_frames:
            message = f"No metric data could be transformed from {metrics_output.path} for run {run_id}."
            logger.error(message)
            raise ValueError(message)

        result = pl.concat(
            output_frames,
            how="diagonal",
        ).with_columns(
            pl.lit(run_id).alias("RUN"),
            pl.lit(metrics_output.workflow_type).alias("WORKFLOW_TYPE"),
            pl.lit(str(metrics_output.workflow_version)).alias("WORKFLOW_VERSION"),
        )

        return self._finalize_frame(result)

    def _standardize_section(
        self,
        section_name: str,
        section: pl.DataFrame,
    ) -> pl.DataFrame:
        """Standardize section columns and metric names."""
        metric_column = section.columns[0]

        if metric_column != self.METRIC_COL:
            section = section.rename({metric_column: (self.METRIC_COL)})

        for required_column in [
            self.LSL_COL,
            self.USL_COL,
        ]:
            if required_column not in section.columns:
                section = section.with_columns(
                    pl.lit(
                        None,
                        dtype=pl.String,
                    ).alias(required_column)
                )

        section = section.with_columns(
            pl.all().cast(
                pl.String,
                strict=False,
            )
        )

        value_columns = [
            column for column in section.columns if column != self.METRIC_COL
        ]

        section = section.with_columns(
            [
                pl.when(
                    pl.col(column).is_null() | pl.col(column).is_in(self.MISSING_VALUES)
                )
                .then(
                    pl.lit(
                        None,
                        dtype=pl.String,
                    )
                )
                .otherwise(pl.col(column).str.strip_chars())
                .alias(column)
                for column in value_columns
            ]
        )

        prefix = self._section_prefix(section_name)

        metric_expression = (
            pl.col(self.METRIC_COL)
            .str.strip_chars()
            .str.replace(
                r"\s+\(.*$",
                "",
            )
            .str.replace_all(
                "%",
                "PCT",
            )
            .str.replace_all(
                r"[^0-9A-Za-z]+",
                "_",
            )
            .str.replace_all(
                r"_+",
                "_",
            )
            .str.strip_chars("_")
        )

        if prefix:
            metric_expression = (
                pl.when(metric_expression.str.starts_with(prefix))
                .then(metric_expression)
                .otherwise(
                    pl.concat_str(
                        [
                            pl.lit(prefix),
                            metric_expression,
                        ]
                    )
                )
            )

        return section.with_columns(metric_expression.alias(self.METRIC_COL)).filter(
            pl.col(self.METRIC_COL).is_not_null() & (pl.col(self.METRIC_COL) != "")
        )

    @staticmethod
    def _section_prefix(
        section_name: str,
    ) -> str:
        """Return the DNA or RNA prefix."""
        upper_name = section_name.upper()

        if "DNA" in upper_name:
            return "DNA_"

        if "RNA" in upper_name:
            return "RNA_"

        return ""

    def _sample_columns(
        self,
        section: pl.DataFrame,
    ) -> list[str]:
        """Return columns representing samples."""
        reserved = {
            self.METRIC_COL,
            self.LSL_COL,
            self.USL_COL,
            self.VALUE_COL,
            "-",
            "",
        }

        return [column for column in section.columns if column not in reserved]

    def _section_to_wide(
        self,
        section: pl.DataFrame,
        value_columns: list[str],
        identifier_mapping: (dict[str, str] | None) = None,
    ) -> pl.DataFrame:
        """Convert one section into a wide table."""
        existing_columns = [
            column for column in value_columns if column in section.columns
        ]

        if not existing_columns:
            return pl.DataFrame()

        mapping = identifier_mapping or {}

        long_frame = (
            section.select(
                [
                    self.METRIC_COL,
                    *existing_columns,
                ]
            )
            .unpivot(
                index=self.METRIC_COL,
                on=existing_columns,
                variable_name="SAMPLE_ID",
                value_name="VALUE",
            )
            .with_columns(pl.col("SAMPLE_ID").replace(mapping).alias("SAMPLE_ID"))
            .filter(pl.col("VALUE").is_not_null())
        )

        if long_frame.is_empty():
            return pl.DataFrame()

        conflicts = (
            long_frame.group_by(
                [
                    "SAMPLE_ID",
                    self.METRIC_COL,
                ]
            )
            .agg(pl.col("VALUE").n_unique().alias("VALUE_COUNT"))
            .filter(pl.col("VALUE_COUNT") > 1)
        )

        if conflicts.height:
            message = (
                f"Conflicting duplicate metric values were detected for run {section}."
            )
            logger.error(message)
            raise ValueError(message)

        return long_frame.pivot(
            index="SAMPLE_ID",
            on=self.METRIC_COL,
            values="VALUE",
            aggregate_function="first",
        )

    def _merge_wide_frames(
        self,
        frames: list[pl.DataFrame],
    ) -> pl.DataFrame:
        """Merge section frames by sample ID."""
        usable_frames = [frame for frame in frames if not frame.is_empty()]

        if not usable_frames:
            return pl.DataFrame()

        combined = pl.concat(
            usable_frames,
            how="diagonal",
        )

        value_columns = [column for column in combined.columns if column != "SAMPLE_ID"]

        for column in value_columns:
            conflicts = (
                combined.group_by("SAMPLE_ID")
                .agg(pl.col(column).drop_nulls().n_unique().alias("VALUE_COUNT"))
                .filter(pl.col("VALUE_COUNT") > 1)
            )

            if conflicts.height:
                message = f"Conflicting values detected for metric {column}."
                logger.error(message)
                raise ValueError(message)

        return combined.group_by(
            "SAMPLE_ID",
            maintain_order=True,
        ).agg(
            [
                pl.col(column).drop_nulls().first().alias(column)
                for column in value_columns
            ]
        )

    def _add_record_type(
        self,
        samples: pl.DataFrame,
    ) -> pl.DataFrame:
        """Classify sample rows as DNA, RNA, or unknown."""
        dna_columns = [
            column for column in samples.columns if column.startswith("DNA_")
        ]
        rna_columns = [
            column for column in samples.columns if column.startswith("RNA_")
        ]

        dna_has_value = (
            pl.any_horizontal([pl.col(column).is_not_null() for column in dna_columns])
            if dna_columns
            else pl.lit(False)
        )

        rna_has_value = (
            pl.any_horizontal([pl.col(column).is_not_null() for column in rna_columns])
            if rna_columns
            else pl.lit(False)
        )

        sample_id_upper = pl.col("SAMPLE_ID").str.to_uppercase()

        dna_id_fallback = sample_id_upper.str.contains(r"^DNA($|_)")

        rna_id_fallback = sample_id_upper.str.contains(r"^RNA($|_)")

        return samples.with_columns(
            pl.when(dna_has_value & ~rna_has_value)
            .then(pl.lit(self.DNA_SAMPLE))
            .when(rna_has_value & ~dna_has_value)
            .then(pl.lit(self.RNA_SAMPLE))
            .when(dna_id_fallback)
            .then(pl.lit(self.DNA_SAMPLE))
            .when(rna_id_fallback)
            .then(pl.lit(self.RNA_SAMPLE))
            .otherwise(pl.lit(self.UNKNOWN_SAMPLE))
            .alias("RECORD_TYPE")
        )

    def _finalize_frame(
        self,
        frame: pl.DataFrame,
    ) -> pl.DataFrame:
        """Order columns and replace null values with NA."""
        metadata = [
            column for column in self.METADATA_COLUMNS if column in frame.columns
        ]

        run_metrics = sorted(
            column
            for column in frame.columns
            if column not in metadata
            and not column.startswith("DNA_")
            and not column.startswith("RNA_")
        )

        dna_metrics = sorted(
            column for column in frame.columns if column.startswith("DNA_")
        )

        rna_metrics = sorted(
            column for column in frame.columns if column.startswith("RNA_")
        )

        return (
            frame.select(metadata + run_metrics + dna_metrics + rna_metrics)
            .with_columns(
                pl.all().cast(
                    pl.String,
                    strict=False,
                )
            )
            .fill_null("NA")
        )

    @staticmethod
    def _combine_runs(
        run_frames: list[pl.DataFrame],
    ) -> pl.DataFrame:
        """Combine all processed workflow frames."""
        if not run_frames:
            logger.error("No run metrics were parsed.")
            raise ValueError("No run metrics were parsed.")

        return pl.concat(
            run_frames,
            how="diagonal",
        ).fill_null("NA")

    def _create_joint_qc(
        self,
        master: pl.DataFrame,
    ) -> pl.DataFrame:
        """Create one joint QC row per run and workflow."""
        sample_rows = master.filter(
            pl.col("RECORD_TYPE").is_in(
                [
                    self.DNA_SAMPLE,
                    self.RNA_SAMPLE,
                    self.UNKNOWN_SAMPLE,
                ]
            )
        )

        aggregations: list[pl.Expr] = []

        for metric in self.JOINT_QC_METRICS:
            if metric in sample_rows.columns:
                aggregations.append(
                    pl.col(metric)
                    .filter(pl.col(metric).is_not_null() & (pl.col(metric) != "NA"))
                    .first()
                    .alias(metric)
                )
            else:
                aggregations.append(
                    pl.lit(
                        None,
                        dtype=pl.String,
                    ).alias(metric)
                )

        return (
            sample_rows.group_by(
                [
                    "RUN_INDEX",
                    "RUN",
                    "WORKFLOW_TYPE",
                    "WORKFLOW_VERSION",
                ],
                maintain_order=True,
            )
            .agg(aggregations)
            .rename(
                {
                    "RUN": "RUN_ID",
                }
            )
            .select(self.JOINT_QC_COLUMNS)
            .fill_null("NA")
        )

    def _add_run_index(
        self,
        frame: pl.DataFrame,
        run_column: str,
    ) -> pl.DataFrame:
        """Add run index where 001 corresponds to the latest run."""

        n_runs = len(self.run_ids)
        width = max(3, len(str(n_runs)))

        mapping = {
            run_id: str(n_runs - index).zfill(width)
            for index, run_id in enumerate(self.run_ids)
        }
        message = (
            f"Assigned run indices to {n_runs} run(s); "
            f"latest run {self.run_ids[-1]} has index 001."
        )
        logger.info(message)
        return frame.with_columns(
            pl.col(run_column).replace(mapping).alias(self.RUN_INDEX)
        )
