"""
metric plots subpackage main module unit tests.
"""

import os
import unittest
from os import path

import polars
import pytest

from tsoppy.metric_plots.main import MetricPlots


# Define path to test data - cannot be absolute due to different paths locally and in CI
test_data_dir = "tests/test_data/metric_plots_main"

config_yaml = "config.yaml"
inpred_nomenclature = path.join(
    test_data_dir,
    "nomenclature.yaml",
)


class TestMetricPlots(unittest.TestCase):
    def test_run(self):
        test_cases = [
            {
                "name": "create master table from dragen workflow",
                "input_directory": path.join(
                    test_data_dir,
                    "dragen_case",
                ),
                "run_ids": [
                    "240809_A02134_0013_BHCGJYDRX5",
                ],
                "workdir": path.join(
                    test_data_dir,
                    "dragen_case",
                    "out",
                ),
                "expected": path.join(
                    test_data_dir,
                    "dragen_case",
                    "master_metrics_table_expected.tsv",
                ),
            },
            {
                "name": "create master table from localapp workflow",
                "input_directory": path.join(
                    test_data_dir,
                    "localapp_case",
                ),
                "run_ids": [
                    "240906_A02134_0019_BHHGKGDRX5",
                ],
                "workdir": path.join(
                    test_data_dir,
                    "localapp_case",
                    "out",
                ),
                "expected": path.join(
                    test_data_dir,
                    "localapp_case",
                    "master_metrics_table_expected.tsv",
                ),
            },
        ]

        for test_case in test_cases:
            with self.subTest(msg=test_case["name"]):
                metric_plots = MetricPlots(
                    config_yaml=config_yaml,
                    inpred_nomenclature=inpred_nomenclature,
                    input_directory=test_case["input_directory"],
                    run_ids=test_case["run_ids"],
                    workdir=test_case["workdir"],
                )

                master, _ = metric_plots.run()

                expected = polars.read_csv(
                    test_case["expected"],
                    separator="\t",
                    infer_schema=False,
                )

                assert master.equals(expected)

                output = path.join(
                    test_case["workdir"],
                    "master_metrics_table.tsv",
                )

                assert path.isfile(output)

                os.remove(output)

                joint_qc = path.join(
                    test_case["workdir"],
                    "joint_sequencing_QC_file.tsv",
                )

                if path.isfile(joint_qc):
                    os.remove(joint_qc)

                if path.isdir(test_case["workdir"]):
                    os.rmdir(test_case["workdir"])

    def test_prepare_plot_frame_last_runs(self):
        master = polars.DataFrame(
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

        metric_plots = MetricPlots.__new__(MetricPlots)

        got = metric_plots.prepare_plot_frame(
            master=master,
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

    def test_prepare_plot_frame_last_runs_more_than_available(self):
        master = polars.DataFrame(
            {
                "SAMPLE_ID": [
                    "S1",
                    "S2",
                    "S3",
                ],
                "RUN": [
                    "RUN1",
                    "RUN2",
                    "RUN3",
                ],
                "WORKFLOW_TYPE": [
                    "dragen",
                    "localapp",
                    "dragen",
                ],
                "WORKFLOW_VERSION": [
                    "2.6.2.4",
                    "ruo-2.2.0.12",
                    "2.6.2.4",
                ],
                "RECORD_TYPE": [
                    "DNA_SAMPLE",
                    "DNA_SAMPLE",
                    "DNA_SAMPLE",
                ],
            }
        )

        metric_plots = MetricPlots.__new__(MetricPlots)

        got = metric_plots.prepare_plot_frame(
            master=master,
            workflow_type="dragen",
            plot_last_runs=10,
        )

        expected = master.filter(polars.col("WORKFLOW_TYPE") == "dragen")

        assert got.equals(expected)

    def test_prepare_plot_frame_explicit_runs(self):
        master = polars.DataFrame(
            {
                "SAMPLE_ID": [
                    "S1",
                    "S2",
                    "S3",
                    "S4",
                ],
                "RUN": [
                    "RUN1",
                    "RUN2",
                    "RUN3",
                    "RUN4",
                ],
                "WORKFLOW_TYPE": [
                    "dragen",
                    "dragen",
                    "localapp",
                    "dragen",
                ],
                "WORKFLOW_VERSION": [
                    "2.6.2.4",
                    "2.6.2.4",
                    "ruo-2.2.0.12",
                    "2.6.2.4",
                ],
                "RECORD_TYPE": [
                    "DNA_SAMPLE",
                    "DNA_SAMPLE",
                    "DNA_SAMPLE",
                    "DNA_SAMPLE",
                ],
            }
        )

        metric_plots = MetricPlots.__new__(MetricPlots)

        got = metric_plots.prepare_plot_frame(
            master=master,
            workflow_type="dragen",
            plot_run_ids=[
                "RUN1",
                "RUN4",
            ],
        )

        expected = master.filter(
            polars.col("RUN").is_in(
                [
                    "RUN1",
                    "RUN4",
                ]
            )
            & (polars.col("WORKFLOW_TYPE") == "dragen")
        )

        assert got.equals(expected)

    def test_prepare_plot_frame_filters_workflow_first(self):
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

        metric_plots = MetricPlots.__new__(MetricPlots)

        got = metric_plots.prepare_plot_frame(
            master=master,
            workflow_type="localapp",
            plot_run_ids=[
                "RUN1",
            ],
        )

        expected = master.filter(polars.col("WORKFLOW_TYPE") == "localapp")
        print("MASTER")
        print(master)
        print(master.schema)

        print("EXPECTED")
        print(expected)
        print(expected.schema)
        assert got.equals(expected)

    def test_prepare_plot_frame_missing_explicit_run(self):
        master = polars.DataFrame(
            {
                "SAMPLE_ID": [
                    "S1",
                ],
                "RUN": [
                    "RUN1",
                ],
                "WORKFLOW_TYPE": [
                    "dragen",
                ],
                "WORKFLOW_VERSION": [
                    "2.6.2.4",
                ],
                "RECORD_TYPE": [
                    "DNA_SAMPLE",
                ],
            }
        )

        metric_plots = MetricPlots.__new__(MetricPlots)

        got = metric_plots.prepare_plot_frame(
            master=master,
            workflow_type="dragen",
            plot_run_ids=[
                "RUN_DOES_NOT_EXIST",
            ],
        )

        assert got.is_empty()
