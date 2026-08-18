"""
Metric plots subpackage main module unit tests.
"""

import os
import unittest
from os import path
import tempfile

import polars

from tsoppy.metric_plots.main import MetricPlots


# Define paths to test data - cannot be absolute due to different paths
# locally and in CI.
test_data_dir = "tests/test_data/metric_plots_main"
config_yaml = "config.yaml"
inpred_nomenclature = "resources/nomenclature.yaml"


class TestMetricPlots(unittest.TestCase):
    @staticmethod
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

    @staticmethod
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

    @staticmethod
    def _metric_plots_without_init():
        """Create MetricPlots without filesystem-dependent initialization."""
        return MetricPlots.__new__(MetricPlots)

    def test_run(self):
        """Create master metrics tables from DRAGEN and LocalApp fixtures."""
        test_cases = [
            {
                "name": "create master table from dragen workflow",
                "input_glob": path.join(
                    test_data_dir,
                    "dragen_case",
                    "dragen",
                    "*",
                ),
                "run_ids": [
                    "240809_A02134_0013_BHCGJYDRX5",
                ],
                "expected": path.join(
                    test_data_dir,
                    "dragen_case",
                    "master_metrics_table_expected.tsv",
                ),
            },
            {
                "name": "create master table from localapp workflow",
                "input_glob": path.join(
                    test_data_dir,
                    "localapp_case",
                    "localapp",
                    "*",
                ),
                "run_ids": ["240906_A02134_0019_BHHGKGDRX5"],
                "expected": path.join(
                    test_data_dir,
                    "localapp_case",
                    "master_metrics_table_expected.tsv",
                ),
            },
        ]

        for test_case in test_cases:
            with self.subTest(msg=test_case["name"]):
                config = path.abspath(config_yaml)
                nomenclature = path.abspath(inpred_nomenclature)
                input_glob = path.abspath(test_case["input_glob"])
                expected_path = path.abspath(test_case["expected"])

                with tempfile.TemporaryDirectory() as tmpdir:
                    current_dir = os.getcwd()

                    try:
                        os.chdir(tmpdir)

                        metric_plots = MetricPlots(
                            config_yaml=config,
                            inpred_nomenclature=nomenclature,
                            input_glob=input_glob,
                            run_ids=test_case["run_ids"],
                        )

                        master, _ = metric_plots.generate_metrics_tables()

                        expected = polars.read_csv(
                            expected_path,
                            separator="\t",
                            infer_schema=False,
                        )

                        assert master.equals(expected)

                        assert path.isfile("master_metrics_table.tsv")

                        assert path.isfile("joint_sequencing_QC_file.tsv")

                    finally:
                        os.chdir(current_dir)

    def test_select_plot_data_last_runs(self):
        master = self._master_frame()
        joint_qc = self._joint_qc_frame()
        metric_plots = self._metric_plots_without_init()

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

    def test_select_plot_data_joint_qc_last_runs(self):
        master = self._master_frame()
        joint_qc = self._joint_qc_frame()
        metric_plots = self._metric_plots_without_init()

        _, got = metric_plots.select_plot_data(
            master=master,
            joint_qc=joint_qc,
            workflow_type="dragen",
            plot_last_runs=1,
        )

        expected = joint_qc.filter(polars.col("RUN_ID") == "RUN5")

        assert got.equals(expected)

    def test_select_plot_data_last_runs_more_than_available(self):
        master = self._master_frame()
        joint_qc = self._joint_qc_frame()
        metric_plots = self._metric_plots_without_init()

        got, _ = metric_plots.select_plot_data(
            master=master,
            joint_qc=joint_qc,
            workflow_type="dragen",
            plot_last_runs=10,
        )

        expected = master.filter(polars.col("WORKFLOW_TYPE") == "dragen")

        assert got.equals(expected)

    def test_select_plot_data_explicit_runs(self):
        master = self._master_frame()
        joint_qc = self._joint_qc_frame()
        metric_plots = self._metric_plots_without_init()

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

    def test_select_plot_data_joint_qc_explicit_runs(self):
        master = self._master_frame()
        joint_qc = self._joint_qc_frame()
        metric_plots = self._metric_plots_without_init()

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

    def test_select_plot_data_filters_workflow_first(self):
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

        metric_plots = self._metric_plots_without_init()

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

    def test_select_plot_data_joint_qc_filters_workflow_first(self):
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

        metric_plots = self._metric_plots_without_init()

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

    def test_select_plot_data_missing_explicit_run(self):
        master = self._master_frame()
        joint_qc = self._joint_qc_frame()
        metric_plots = self._metric_plots_without_init()

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

    def test_select_plot_data_filters_both_outputs_consistently(self):
        master = self._master_frame()
        joint_qc = self._joint_qc_frame()
        metric_plots = self._metric_plots_without_init()

        plot_frame, plot_joint_qc = metric_plots.select_plot_data(
            master=master,
            joint_qc=joint_qc,
            workflow_type="dragen",
            plot_last_runs=1,
        )

        assert plot_frame["RUN"].to_list() == ["RUN5"]

        assert plot_joint_qc["RUN_ID"].to_list() == ["RUN5"]
