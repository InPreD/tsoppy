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


def test_select_plot_data_last_runs():
    master = _master_frame()
    joint_qc = _joint_qc_frame()
    metric_plots = _metric_plots_without_init()

    got, _ = metric_plots.select_plot_data(
        master=master,
        joint_qc=joint_qc,
        workflow_type="dragen",
        plot_last_runs=2,
    )

    expected = master.filter(
        polars.col("RUN").is_in(
            [
                "RUN3",
                "RUN5",
            ]
        )
    )

    assert got.equals(expected)


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


def test_select_plot_data_defaults_to_last_ten_runs():
    """Select the last ten workflow runs when no plot selector is provided."""

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

    metric_plots = _metric_plots_without_init()

    got, got_joint_qc = metric_plots.select_plot_data(
        master=master,
        joint_qc=joint_qc,
        workflow_type="dragen",
    )

    expected_runs = [f"RUN{i:02d}" for i in range(3, 13)]

    assert got["RUN"].to_list() == expected_runs
    assert got_joint_qc["RUN_ID"].to_list() == expected_runs


def test_select_plot_data_joint_qc_last_runs():
    master = _master_frame()
    joint_qc = _joint_qc_frame()
    metric_plots = _metric_plots_without_init()

    _, got = metric_plots.select_plot_data(
        master=master,
        joint_qc=joint_qc,
        workflow_type="dragen",
        plot_last_runs=1,
    )

    expected = joint_qc.filter(polars.col("RUN_ID") == "RUN5")

    assert got.equals(expected)


def test_select_plot_data_last_runs_more_than_available():
    master = _master_frame()
    joint_qc = _joint_qc_frame()
    metric_plots = _metric_plots_without_init()

    got, _ = metric_plots.select_plot_data(
        master=master,
        joint_qc=joint_qc,
        workflow_type="dragen",
        plot_last_runs=10,
    )

    expected = master.filter(polars.col("WORKFLOW_TYPE") == "dragen")

    assert got.equals(expected)


def test_select_plot_data_explicit_runs():
    master = _master_frame()
    joint_qc = _joint_qc_frame()
    metric_plots = _metric_plots_without_init()

    got, _ = metric_plots.select_plot_data(
        master=master,
        joint_qc=joint_qc,
        workflow_type="dragen",
        plot_run_ids=[
            "RUN1",
            "RUN5",
        ],
    )

    expected = master.filter(
        polars.col("RUN").is_in(
            [
                "RUN1",
                "RUN5",
            ]
        )
        & (polars.col("WORKFLOW_TYPE") == "dragen")
    )

    assert got.equals(expected)


def test_select_plot_data_joint_qc_explicit_runs():
    master = _master_frame()
    joint_qc = _joint_qc_frame()
    metric_plots = _metric_plots_without_init()

    _, got = metric_plots.select_plot_data(
        master=master,
        joint_qc=joint_qc,
        workflow_type="dragen",
        plot_run_ids=[
            "RUN1",
            "RUN5",
        ],
    )

    expected = joint_qc.filter(
        polars.col("RUN_ID").is_in(
            [
                "RUN1",
                "RUN5",
            ]
        )
        & (polars.col("WORKFLOW_TYPE") == "dragen")
    )

    assert got.equals(expected)


def test_select_plot_data_filters_workflow_first():
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

    got, _ = metric_plots.select_plot_data(
        master=master,
        joint_qc=joint_qc,
        workflow_type="localapp",
        plot_run_ids=[
            "RUN1",
        ],
    )

    expected = master.filter(polars.col("WORKFLOW_TYPE") == "localapp")

    assert got.equals(expected)


def test_select_plot_data_joint_qc_filters_workflow_first():
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

    _, got = metric_plots.select_plot_data(
        master=master,
        joint_qc=joint_qc,
        workflow_type="localapp",
        plot_run_ids=[
            "RUN1",
        ],
    )

    expected = joint_qc.filter(polars.col("WORKFLOW_TYPE") == "localapp")

    assert got.equals(expected)


def test_select_plot_data_missing_explicit_run():
    master = _master_frame()
    joint_qc = _joint_qc_frame()
    metric_plots = _metric_plots_without_init()

    got, got_joint_qc = metric_plots.select_plot_data(
        master=master,
        joint_qc=joint_qc,
        workflow_type="dragen",
        plot_run_ids=[
            "RUN_DOES_NOT_EXIST",
        ],
    )

    assert got.is_empty()
    assert got_joint_qc.is_empty()


def test_select_plot_data_filters_both_outputs_consistently():
    master = _master_frame()
    joint_qc = _joint_qc_frame()
    metric_plots = _metric_plots_without_init()

    plot_frame, plot_joint_qc = metric_plots.select_plot_data(
        master=master,
        joint_qc=joint_qc,
        workflow_type="dragen",
        plot_last_runs=1,
    )

    assert plot_frame["RUN"].to_list() == ["RUN5"]

    assert plot_joint_qc["RUN_ID"].to_list() == ["RUN5"]


def test_add_record_type_uses_samplesheet_sample_type():
    """RECORD_TYPE comes from the sample sheet, not metric content or SAMPLE_ID text."""
    samples = polars.DataFrame(
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
    )

    samplesheet = polars.DataFrame(
        {
            "Sample_ID": [
                "RNA_LOOKING_ID",
                "DNA_LOOKING_ID",
            ],
            "Pair_ID": [
                "RNA_LOOKING_ID",
                "DNA_LOOKING_ID",
            ],
            "Sample_Type": [
                "DNA",
                "RNA",
            ],
        }
    )

    metric_plots = _metric_plots_without_init()

    got = metric_plots._add_record_type(samples, samplesheet)

    assert got["RECORD_TYPE"].to_list() == [
        "DNA_SAMPLE",
        "RNA_SAMPLE",
    ]


def test_add_record_type_prefers_pair_id_over_sample_id():
    """The metrics output keys samples by Pair_ID, so the lookup joins on it."""
    samples = polars.DataFrame(
        {
            "SAMPLE_ID": ["PAIR01"],
        }
    )

    samplesheet = polars.DataFrame(
        {
            "Sample_ID": ["SAMPLE01"],
            "Pair_ID": ["PAIR01"],
            "Sample_Type": ["RNA"],
        }
    )

    metric_plots = _metric_plots_without_init()

    got = metric_plots._add_record_type(samples, samplesheet)

    assert got["RECORD_TYPE"].to_list() == ["RNA_SAMPLE"]


def test_add_record_type_falls_back_to_sample_id_without_pair_id_column():
    samples = polars.DataFrame(
        {
            "SAMPLE_ID": ["SAMPLE01"],
        }
    )

    samplesheet = polars.DataFrame(
        {
            "Sample_ID": ["SAMPLE01"],
            "Sample_Type": ["DNA"],
        }
    )

    metric_plots = _metric_plots_without_init()

    got = metric_plots._add_record_type(samples, samplesheet)

    assert got["RECORD_TYPE"].to_list() == ["DNA_SAMPLE"]


def test_add_record_type_unmatched_sample_falls_back_to_unknown():
    samples = polars.DataFrame(
        {
            "SAMPLE_ID": ["NOT_IN_SAMPLESHEET"],
        }
    )

    samplesheet = polars.DataFrame(
        {
            "Sample_ID": ["OTHER_SAMPLE"],
            "Pair_ID": ["OTHER_SAMPLE"],
            "Sample_Type": ["DNA"],
        }
    )

    metric_plots = _metric_plots_without_init()

    got = metric_plots._add_record_type(samples, samplesheet)

    assert got["RECORD_TYPE"].to_list() == ["SAMPLE"]


def test_add_record_type_missing_sample_type_column_falls_back_to_unknown():
    samples = polars.DataFrame(
        {
            "SAMPLE_ID": ["SAMPLE01"],
        }
    )

    samplesheet = polars.DataFrame(
        {
            "Sample_ID": ["SAMPLE01"],
            "Pair_ID": ["SAMPLE01"],
        }
    )

    metric_plots = _metric_plots_without_init()

    got = metric_plots._add_record_type(samples, samplesheet)

    assert got["RECORD_TYPE"].to_list() == ["SAMPLE"]


def test_add_record_type_ambiguous_samplesheet_uses_fallback():
    """A Pair_ID shared by a DNA and an RNA sample cannot be classified unambiguously."""
    samples = polars.DataFrame(
        {
            "SAMPLE_ID": [
                "SHARED_PAIR",
                "SOLO_PAIR",
            ],
        }
    )

    samplesheet = polars.DataFrame(
        {
            "Sample_ID": [
                "Patient01_D",
                "Patient01_R",
                "Patient02",
            ],
            "Pair_ID": [
                "SHARED_PAIR",
                "SHARED_PAIR",
                "SOLO_PAIR",
            ],
            "Sample_Type": [
                "DNA",
                "RNA",
                "RNA",
            ],
        }
    )

    metric_plots = _metric_plots_without_init()

    got = metric_plots._add_record_type(samples, samplesheet)

    assert got["RECORD_TYPE"].to_list() == [
        "SAMPLE",
        "RNA_SAMPLE",
    ]


def test_add_record_type_sample_type_is_case_and_whitespace_insensitive():
    samples = polars.DataFrame(
        {
            "SAMPLE_ID": ["SAMPLE01"],
        }
    )

    samplesheet = polars.DataFrame(
        {
            "Sample_ID": ["SAMPLE01"],
            "Pair_ID": ["SAMPLE01"],
            "Sample_Type": [" dna "],
        }
    )

    metric_plots = _metric_plots_without_init()

    got = metric_plots._add_record_type(samples, samplesheet)

    assert got["RECORD_TYPE"].to_list() == ["DNA_SAMPLE"]
