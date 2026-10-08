"""Metric plots subpackage main module unit tests."""

import os
import tempfile
from os import path

import polars
import pytest

from tsoppy.metric_plots.main import MetricPlots

# Define paths to test data - cannot be absolute due to different paths
# locally and in CI.
test_data_dir = "tests/test_data/metric_plots_main"
config_yaml = "config.yaml"
inpred_nomenclature = "tests/test_data/metric_plots_main/nomenclature.yaml"


def _master_frame():
    """Create a synthetic master metrics dataframe for plotting tests."""
    return polars.DataFrame(
        {
            "SAMPLE_ID": [
                "S1",
                "S2",
                "S3",
                "S4",
                "S5",
            ],
            "RUN": [
                "RUN1",
                "RUN2",
                "RUN3",
                "RUN4",
                "RUN5",
            ],
            "WORKFLOW_TYPE": [
                "dragen",
                "localapp",
                "dragen",
                "localapp",
                "dragen",
            ],
            "WORKFLOW_VERSION": [
                "2.6.2.4",
                "ruo-2.2.0.12",
                "2.6.2.4",
                "ruo-2.2.0.12",
                "2.6.2.4",
            ],
            "RECORD_TYPE": [
                "DNA_SAMPLE",
                "DNA_SAMPLE",
                "DNA_SAMPLE",
                "DNA_SAMPLE",
                "DNA_SAMPLE",
            ],
        }
    )


def _joint_qc_frame():
    """Create a synthetic joint QC dataframe for plotting tests."""
    return polars.DataFrame(
        {
            "RUN_ID": [
                "RUN1",
                "RUN2",
                "RUN3",
                "RUN4",
                "RUN5",
            ],
            "WORKFLOW_TYPE": [
                "dragen",
                "localapp",
                "dragen",
                "localapp",
                "dragen",
            ],
            "WORKFLOW_VERSION": [
                "2.6.2.4",
                "ruo-2.2.0.12",
                "2.6.2.4",
                "ruo-2.2.0.12",
                "2.6.2.4",
            ],
        }
    )


def _twelve_run_frames():
    """Create 12 same-workflow runs to exercise the default last-10-runs cutoff."""
    run_ids = [f"RUN{i:02d}" for i in range(1, 13)]

    master = polars.DataFrame(
        {
            "SAMPLE_ID": [f"S{i:02d}" for i in range(1, 13)],
            "RUN": run_ids,
            "WORKFLOW_TYPE": ["dragen"] * 12,
            "WORKFLOW_VERSION": ["2.6.2.4"] * 12,
            "RECORD_TYPE": ["DNA_SAMPLE"] * 12,
        }
    )

    joint_qc = polars.DataFrame(
        {
            "RUN_ID": run_ids,
            "WORKFLOW_TYPE": ["dragen"] * 12,
            "WORKFLOW_VERSION": ["2.6.2.4"] * 12,
        }
    )

    return master, joint_qc


def _metric_plots_without_init():
    """Create MetricPlots without filesystem-dependent initialization."""
    return MetricPlots.__new__(MetricPlots)


@pytest.mark.parametrize(
    ("input_glob", "run_ids", "expected"),
    [
        (
            # create master table from dragen workflow
            path.join(
                test_data_dir,
                "dragen_case",
                "dragen",
                "*",
            ),
            [
                "240809_A02134_0013_BHCGJYDRX5",
            ],
            path.join(
                test_data_dir,
                "dragen_case",
                "master_metrics_table_expected.tsv",
            ),
        ),
        (
            # create master table from localapp workflow
            path.join(
                test_data_dir,
                "localapp_case",
                "localapp",
                "*",
            ),
            ["240906_A02134_0019_BHHGKGDRX5"],
            path.join(
                test_data_dir,
                "localapp_case",
                "master_metrics_table_expected.tsv",
            ),
        ),
    ],
)
def test_run(input_glob, run_ids, expected):
    """Create master metrics tables from DRAGEN and LocalApp fixtures."""
    config = path.abspath(config_yaml)
    nomenclature = path.abspath(inpred_nomenclature)
    input_glob = path.abspath(input_glob)
    expected_path = path.abspath(expected)

    with tempfile.TemporaryDirectory() as tmpdir:
        current_dir = os.getcwd()

        try:
            os.chdir(tmpdir)

            metric_plots = MetricPlots(
                config_yaml=config,
                inpred_nomenclature=nomenclature,
                input_glob=input_glob,
                run_ids=run_ids,
            )

            master, _ = metric_plots.generate_metrics_tables()

            want = polars.read_csv(
                expected_path,
                separator="\t",
                infer_schema=False,
            )

            assert master.equals(want)

            assert path.isfile("master_metrics_table.tsv")

            assert path.isfile("joint_sequencing_QC_file.tsv")

        finally:
            os.chdir(current_dir)


@pytest.mark.parametrize(
    ("master", "joint_qc", "plot_last_runs", "plot_run_ids", "want_runs"),
    [
        (
            # last_runs selects the N most recent runs for the workflow
            _master_frame(),
            _joint_qc_frame(),
            2,
            None,
            ["RUN3", "RUN5"],
        ),
        (
            # last_runs beyond what's available clips to all available runs
            _master_frame(),
            _joint_qc_frame(),
            10,
            None,
            ["RUN1", "RUN3", "RUN5"],
        ),
        (
            # explicit run IDs are selected regardless of recency
            _master_frame(),
            _joint_qc_frame(),
            None,
            ["RUN1", "RUN5"],
            ["RUN1", "RUN5"],
        ),
        (
            # an explicit run ID matching nothing selects no rows
            _master_frame(),
            _joint_qc_frame(),
            None,
            ["RUN_DOES_NOT_EXIST"],
            [],
        ),
        (
            # no selector defaults to the last DEFAULT_PLOT_LAST_RUNS (10) runs
            *_twelve_run_frames(),
            None,
            None,
            [f"RUN{i:02d}" for i in range(3, 13)],
        ),
    ],
)
def test_select_plot_data(
    master,
    joint_qc,
    plot_last_runs,
    plot_run_ids,
    want_runs,
):
    """Master and joint-QC rows are both filtered to the same selected runs."""
    metric_plots = _metric_plots_without_init()

    got, got_joint_qc = metric_plots.select_plot_data(
        master=master,
        joint_qc=joint_qc,
        workflow_type="dragen",
        plot_last_runs=plot_last_runs,
        plot_run_ids=plot_run_ids,
    )

    assert got["RUN"].to_list() == want_runs
    assert got_joint_qc["RUN_ID"].to_list() == want_runs


def test_no_run_selector_uses_all_runs_from_input_glob():
    """Use all glob-matched run IDs when no run selector is provided."""

    with tempfile.TemporaryDirectory() as tmpdir:
        for run_id in [
            "RUN003",
            "RUN001",
            "RUN002",
        ]:
            os.mkdir(path.join(tmpdir, run_id))

        metric_plots = MetricPlots(
            config_yaml=config_yaml,
            inpred_nomenclature=inpred_nomenclature,
            input_glob=path.join(tmpdir, "*"),
        )

        assert metric_plots.run_ids == [
            "RUN001",
            "RUN002",
            "RUN003",
        ]


def test_select_plot_data_filters_workflow_first():
    """Shared run IDs across workflows do not leak between master and joint QC."""
    master = polars.DataFrame(
        {
            "SAMPLE_ID": [
                "DRAGEN_SAMPLE",
                "LOCALAPP_SAMPLE",
            ],
            "RUN": [
                "RUN1",
                "RUN1",
            ],
            "WORKFLOW_TYPE": [
                "dragen",
                "localapp",
            ],
            "WORKFLOW_VERSION": [
                "2.6.2.4",
                "ruo-2.2.0.12",
            ],
            "RECORD_TYPE": [
                "DNA_SAMPLE",
                "DNA_SAMPLE",
            ],
        }
    )

    joint_qc = polars.DataFrame(
        {
            "RUN_ID": [
                "RUN1",
                "RUN1",
            ],
            "WORKFLOW_TYPE": [
                "dragen",
                "localapp",
            ],
            "WORKFLOW_VERSION": [
                "2.6.2.4",
                "ruo-2.2.0.12",
            ],
        }
    )

    metric_plots = _metric_plots_without_init()

    got, got_joint_qc = metric_plots.select_plot_data(
        master=master,
        joint_qc=joint_qc,
        workflow_type="localapp",
        plot_run_ids=[
            "RUN1",
        ],
    )

    assert got.equals(master.filter(polars.col("WORKFLOW_TYPE") == "localapp"))

    assert got_joint_qc.equals(
        joint_qc.filter(polars.col("WORKFLOW_TYPE") == "localapp")
    )


@pytest.mark.parametrize(
    ("samples", "samplesheet", "want"),
    [
        (
            # RECORD_TYPE comes from the sample sheet, not metric content or SAMPLE_ID text
            polars.DataFrame(
                {
                    "SAMPLE_ID": [
                        "RNA_LOOKING_ID",
                        "DNA_LOOKING_ID",
                    ],
                    "DNA_METRIC": [
                        "10",
                        None,
                    ],
                    "RNA_METRIC": [
                        None,
                        "20",
                    ],
                }
            ),
            polars.DataFrame(
                {
                    "Sample_ID": [
                        "RNA_LOOKING_ID",
                        "DNA_LOOKING_ID",
                    ],
                    "Sample_Type": [
                        "DNA",
                        "RNA",
                    ],
                }
            ),
            [
                "DNA_SAMPLE",
                "RNA_SAMPLE",
            ],
        ),
        (
            # a sample absent from the sample sheet falls back to unknown
            polars.DataFrame(
                {
                    "SAMPLE_ID": ["NOT_IN_SAMPLESHEET"],
                }
            ),
            polars.DataFrame(
                {
                    "Sample_ID": ["OTHER_SAMPLE"],
                    "Sample_Type": ["DNA"],
                }
            ),
            ["SAMPLE"],
        ),
        (
            # a sample sheet without a Sample_Type column falls back to unknown
            polars.DataFrame(
                {
                    "SAMPLE_ID": ["SAMPLE01"],
                }
            ),
            polars.DataFrame(
                {
                    "Sample_ID": ["SAMPLE01"],
                }
            ),
            ["SAMPLE"],
        ),
        (
            # a Sample_ID mapped to more than one Sample_Type cannot be classified unambiguously
            polars.DataFrame(
                {
                    "SAMPLE_ID": [
                        "DUPLICATE_ID",
                        "SOLO_SAMPLE",
                    ],
                }
            ),
            polars.DataFrame(
                {
                    "Sample_ID": [
                        "DUPLICATE_ID",
                        "DUPLICATE_ID",
                        "SOLO_SAMPLE",
                    ],
                    "Sample_Type": [
                        "DNA",
                        "RNA",
                        "RNA",
                    ],
                }
            ),
            [
                "SAMPLE",
                "RNA_SAMPLE",
            ],
        ),
        (
            # Sample_Type matching is case- and whitespace-insensitive
            polars.DataFrame(
                {
                    "SAMPLE_ID": ["SAMPLE01"],
                }
            ),
            polars.DataFrame(
                {
                    "Sample_ID": ["SAMPLE01"],
                    "Sample_Type": [" dna "],
                }
            ),
            ["DNA_SAMPLE"],
        ),
    ],
)
def test_add_record_type(samples, samplesheet, want):
    metric_plots = _metric_plots_without_init()

    got = metric_plots._add_record_type(samples, samplesheet)

    assert got["RECORD_TYPE"].to_list() == want
