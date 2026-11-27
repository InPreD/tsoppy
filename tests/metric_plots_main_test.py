"""Metric plots subpackage main module unit tests."""

import os
from os import path

import polars as pl
import pytest

from tsoppy.metric_plots.main import MetricPlots

# Define paths to test data - cannot be absolute due to different paths
# locally and in CI.
test_data_dir = "tests/test_data/metric_plots_main"
config_yaml = "config.yaml"
inpred_nomenclature = path.join(test_data_dir, "nomenclature.yaml")

# Test data representing master metrics table
test_sample_id_column = [
    "LSL_Guideline",
    "USL_Guideline",
    "DNA_TEST",
    "RNA_TEST",
]

test_record_type_column = [
    "LOWER_THRESHOLD",
    "UPPER_THRESHOLD",
    "DNA_SAMPLE",
    "RNA_SAMPLE",
]

test_completed_all_steps_column = [
    "NA",
    "NA",
    "TRUE",
    "TRUE",
]

test_dna_contamination_score_column = [
    "0",
    "1457",
    "300",
    "NA",
]

test_dna_mad_column = [
    "0.000",
    "0.134",
    "0.066",
    "NA",
]

test_dna_mean_family_size_column = [
    "NA",
    "NA",
    "1.2",
    "NA",
]

test_dna_mean_target_coverage_column = [
    "NA",
    "NA",
    "929.2",
    "NA",
]

test_dna_median_bin_count_cnv_target_column = [
    "1.0",
    "NA",
    "9.3",
    "NA",
]

test_dna_median_exon_coverage_column = [
    "150",
    "NA",
    "1024",
    "NA",
]

test_dna_median_insert_size_column = [
    "70",
    "NA",
    "148",
    "NA",
]

test_dna_median_target_coverage_column = [
    "NA",
    "NA",
    "947.0",
    "NA",
]

test_dna_pct_aligned_reads_column = [
    "NA",
    "NA",
    "99.3",
    "NA",
]

test_dna_pct_chimeric_reads_column = [
    "NA",
    "8.00",
    "0.11",
    "NA",
]

test_dna_pct_contamination_est_column = [
    "NA",
    "NA",
    "17.0",
    "NA",
]

test_dna_pct_exon_100x_column = [
    "NA",
    "NA",
    "99.8",
    "NA",
]

test_dna_pct_exon_50x_column = [
    "90.0",
    "NA",
    "99.8",
    "NA",
]

test_dna_pct_read_enrichment_column = [
    "NA",
    "NA",
    "63.2",
    "NA",
]

test_dna_pct_target_0_4x_mean_column = [
    "NA",
    "NA",
    "92.7",
    "NA",
]

test_dna_pct_target_100x_column = [
    "NA",
    "NA",
    "99.4",
    "NA",
]

test_dna_pct_target_250x_column = [
    "NA",
    "NA",
    "97.0",
    "NA",
]

test_dna_pct_usable_umi_reads_column = [
    "NA",
    "NA",
    "99.7",
    "NA",
]

test_dna_total_pf_reads_column = [
    "NA",
    "NA",
    "95960542",
    "NA",
]

test_dna_usable_msi_sites_column = [
    "40",
    "NA",
    "117",
    "NA",
]

test_rna_median_cv_gene_500x_column = [
    "0.00",
    "0.93",
    "NA",
    "0.48",
]

test_rna_median_insert_size_column = [
    "80",
    "NA",
    "NA",
    "150",
]

test_rna_pct_chimeric_reads_column = ["NA"] * 3 + ["0.46"]

test_rna_pct_on_target_reads_column = ["NA"] * 3 + ["71.4"]

test_rna_scaled_median_gene_coverage_column = ["NA"] * 3 + ["2682.5"]

test_rna_total_on_target_reads_column = [
    "2500000",
    "NA",
    "NA",
    "11889549",
]

test_rna_total_pf_reads_column = ["NA"] * 3 + ["16645676"]

test_master_metrics_table_dragen = pl.DataFrame(
    {
        "RUN_INDEX": ["001"] * 4,
        "SAMPLE_ID": test_sample_id_column,
        "RUN": ["200101_A01234_0001_ABCDE12345"] * 4,
        "WORKFLOW_TYPE": ["dragen"] * 4,
        "WORKFLOW_VERSION": ["2.6.2.4"] * 4,
        "RECORD_TYPE": test_record_type_column,
        "COMPLETED_ALL_STEPS": test_completed_all_steps_column,
        "DNA_ALLELE_DOSAGE_RATIO": [
            "NA",
            "NA",
            "50.000",
            "NA",
        ],
        "DNA_CONTAMINATION_SCORE": test_dna_contamination_score_column,
        "DNA_EXCESSIVE_TF": [
            "NA",
            "0",
            "0",
            "NA",
        ],
        "DNA_GENE_SCALED_MAD": test_dna_mad_column,
        "DNA_MEAN_FAMILY_SIZE": test_dna_mean_family_size_column,
        "DNA_MEAN_TARGET_COVERAGE": test_dna_mean_target_coverage_column,
        "DNA_MEDIAN_BIN_COUNT_CNV_TARGET": test_dna_median_bin_count_cnv_target_column,
        "DNA_MEDIAN_EXON_COVERAGE": test_dna_median_exon_coverage_column,
        "DNA_MEDIAN_INSERT_SIZE": test_dna_median_insert_size_column,
        "DNA_MEDIAN_TARGET_COVERAGE": test_dna_median_target_coverage_column,
        "DNA_MEDIAN_TARGET_HRD_COVERAGE": [
            "NA",
            "NA",
            "284.0",
            "NA",
        ],
        "DNA_PCT_ALIGNED_READS": test_dna_pct_aligned_reads_column,
        "DNA_PCT_CHIMERIC_READS": test_dna_pct_chimeric_reads_column,
        "DNA_PCT_CONTAMINATION_EST": test_dna_pct_contamination_est_column,
        "DNA_PCT_EXON_100X": test_dna_pct_exon_100x_column,
        "DNA_PCT_EXON_50X": test_dna_pct_exon_50x_column,
        "DNA_PCT_Q30_BASES": [
            "NA",
            "NA",
            "90.95",
            "NA",
        ],
        "DNA_PCT_READ_ENRICHMENT": test_dna_pct_read_enrichment_column,
        "DNA_PCT_SOFT_CLIPPED_BASES": [
            "NA",
            "NA",
            "2.16",
            "NA",
        ],
        "DNA_PCT_TARGET_0_4X_MEAN": test_dna_pct_target_0_4x_mean_column,
        "DNA_PCT_TARGET_100X": test_dna_pct_target_100x_column,
        "DNA_PCT_TARGET_250X": test_dna_pct_target_250x_column,
        "DNA_PCT_TARGET_50X": [
            "NA",
            "NA",
            "99.7",
            "NA",
        ],
        "DNA_PCT_TARGET_HRD_50X": [
            "50.0",
            "NA",
            "99.6",
            "NA",
        ],
        "DNA_PCT_USABLE_UMI_READS": test_dna_pct_usable_umi_reads_column,
        "DNA_TOTAL_PF_READS": test_dna_total_pf_reads_column,
        "DNA_USABLE_MSI_SITES": test_dna_usable_msi_sites_column,
        "RNA_GENE_ABOVE_MEDIAN_CUTOFF": ["NA"] * 3 + ["42"],
        "RNA_GENE_MEDIAN_COVERAGE": ["NA"] * 3 + ["966.24"],
        "RNA_MEDIAN_CV_GENE_500X": test_rna_median_cv_gene_500x_column,
        "RNA_MEDIAN_INSERT_SIZE": test_rna_median_insert_size_column,
        "RNA_PCT_CHIMERIC_READS": test_rna_pct_chimeric_reads_column,
        "RNA_PCT_ON_TARGET_READS": test_rna_pct_on_target_reads_column,
        "RNA_PCT_Q30_BASES": ["NA"] * 3 + ["90.60"],
        "RNA_PCT_SOFT_CLIPPED_BASES": ["NA"] * 3 + ["2.18"],
        "RNA_SCALED_MEDIAN_GENE_COVERAGE": test_rna_scaled_median_gene_coverage_column,
        "RNA_TOTAL_ON_TARGET_READS": test_rna_total_on_target_reads_column,
        "RNA_TOTAL_PF_READS": test_rna_total_pf_reads_column,
    }
)

test_pct_q30_column = [
    "80.0",
    "NA",
    "95",
    "95",
]

test_master_metrics_table_localapp = pl.DataFrame(
    {
        "RUN_INDEX": ["001"] * 4,
        "SAMPLE_ID": test_sample_id_column,
        "RUN": ["200101_A01234_0002_ABCDE12345"] * 4,
        "WORKFLOW_TYPE": ["localapp"] * 4,
        "WORKFLOW_VERSION": ["ruo-2.2.0.12"] * 4,
        "RECORD_TYPE": test_record_type_column,
        "COMPLETED_ALL_STEPS": test_completed_all_steps_column,
        "PCT_PF_READS": [
            "55.0",
            "NA",
            "89.1",
            "89.1",
        ],
        "PCT_Q30_R1": test_pct_q30_column,
        "PCT_Q30_R2": test_pct_q30_column,
        "DNA_CONTAMINATION_P_VALUE": [
            "0.000",
            "0.049",
            "1.000",
            "NA",
        ],
        "DNA_CONTAMINATION_SCORE": test_dna_contamination_score_column,
        "DNA_COVERAGE_MAD": test_dna_mad_column,
        "DNA_MEAN_FAMILY_SIZE": test_dna_mean_family_size_column,
        "DNA_MEAN_TARGET_COVERAGE": test_dna_mean_target_coverage_column,
        "DNA_MEDIAN_BIN_COUNT_CNV_TARGET": test_dna_median_bin_count_cnv_target_column,
        "DNA_MEDIAN_EXON_COVERAGE": test_dna_median_exon_coverage_column,
        "DNA_MEDIAN_INSERT_SIZE": test_dna_median_insert_size_column,
        "DNA_MEDIAN_TARGET_COVERAGE": test_dna_median_target_coverage_column,
        "DNA_PCT_ALIGNED_READS": test_dna_pct_aligned_reads_column,
        "DNA_PCT_CHIMERIC_READS": test_dna_pct_chimeric_reads_column,
        "DNA_PCT_CONTAMINATION_EST": test_dna_pct_contamination_est_column,
        "DNA_PCT_EXON_100X": test_dna_pct_exon_100x_column,
        "DNA_PCT_EXON_50X": test_dna_pct_exon_50x_column,
        "DNA_PCT_PF_UQ_READS": [
            "NA",
            "NA",
            "100.0",
            "NA",
        ],
        "DNA_PCT_READ_ENRICHMENT": test_dna_pct_read_enrichment_column,
        "DNA_PCT_TARGET_0_4X_MEAN": test_dna_pct_target_0_4x_mean_column,
        "DNA_PCT_TARGET_100X": test_dna_pct_target_100x_column,
        "DNA_PCT_TARGET_250X": test_dna_pct_target_250x_column,
        "DNA_PCT_USABLE_UMI_READS": test_dna_pct_usable_umi_reads_column,
        "DNA_TOTAL_PF_READS": test_dna_total_pf_reads_column,
        "DNA_USABLE_MSI_SITES": test_dna_usable_msi_sites_column,
        "RNA_MEDIAN_CV_GENE_500X": test_rna_median_cv_gene_500x_column,
        "RNA_MEDIAN_INSERT_SIZE": test_rna_median_insert_size_column,
        "RNA_PCT_CHIMERIC_READS": test_rna_pct_chimeric_reads_column,
        "RNA_PCT_ON_TARGET_READS": test_rna_pct_on_target_reads_column,
        "RNA_SCALED_MEDIAN_GENE_COVERAGE": test_rna_scaled_median_gene_coverage_column,
        "RNA_TOTAL_ON_TARGET_READS": test_rna_total_on_target_reads_column,
        "RNA_TOTAL_PF_READS": test_rna_total_pf_reads_column,
    }
)

# Test data representing master frame and joint qc frame
test_workflow_type_column = ["dragen", "localapp"] * 6

test_workflow_version_column = ["2.6.2.4", "ruo-2.2.0.12"] * 6

test_master_frame = pl.DataFrame(
    {
        "SAMPLE_ID": [f"S{var}" for var in range(1, 13)],
        "RUN": [f"RUN{var}" for var in range(1, 13)],
        "WORKFLOW_TYPE": test_workflow_type_column,
        "WORKFLOW_VERSION": test_workflow_version_column,
        "RECORD_TYPE": ["DNA_SAMPLE"] * 12,
    }
)

test_joint_qc_frame = pl.DataFrame(
    {
        "RUN_ID": [f"RUN{var}" for var in range(1, 13)],
        "WORKFLOW_TYPE": test_workflow_type_column,
        "WORKFLOW_VERSION": test_workflow_version_column,
    }
)


def test_metricplots():
    """Use all glob-matched run IDs when no run selector is provided."""
    metric_plots = MetricPlots(
        config_yaml=config_yaml,
        inpred_nomenclature=inpred_nomenclature,
        input_glob=path.join(
            test_data_dir,
            "metricplots",
            "*",
        ),
    )
    assert metric_plots.run_ids == [
        "RUN001",
        "RUN002",
        "RUN003",
    ]


@pytest.mark.parametrize(
    ("inputs, want"),
    [
        (
            # create master table from dragen workflow
            (
                path.join(
                    test_data_dir,
                    "dragen",
                    "*",
                ),
                ["200101_A01234_0001_ABCDE12345"],
            ),
            test_master_metrics_table_dragen,
        ),
        (
            # create master table from localapp workflow
            (
                path.join(
                    test_data_dir,
                    "localapp",
                    "*",
                ),
                ["200101_A01234_0002_ABCDE12345"],
            ),
            test_master_metrics_table_localapp,
        ),
    ],
)
def test_metricplots_generate_metrics_tables(inputs, want):
    """Create master metrics tables from DRAGEN and LocalApp fixtures."""
    metric_plots = MetricPlots(
        config_yaml=config_yaml,
        inpred_nomenclature=inpred_nomenclature,
        input_glob=inputs[0],
        run_ids=inputs[1],
    )
    got, _ = metric_plots.generate_metrics_tables()
    assert got.equals(want)
    for got_file in ["master_metrics_table.tsv", "joint_sequencing_QC_file.tsv"]:
        assert path.isfile(got_file)
        os.remove(got_file)


@pytest.mark.parametrize(
    ("inputs, wants"),
    [
        (
            # last_runs selects the N most recent runs for the workflow
            (
                test_master_frame,
                test_joint_qc_frame,
                "dragen",
                2,
                None,
            ),
            (
                pl.concat(
                    [test_master_frame.slice(8, 1), test_master_frame.slice(10, 1)]
                ),
                pl.concat(
                    [test_joint_qc_frame.slice(8, 1), test_joint_qc_frame.slice(10, 1)]
                ),
            ),
        ),
        (
            # explicit run IDs are selected regardless of recency
            (
                test_master_frame,
                test_joint_qc_frame,
                "dragen",
                None,
                ["RUN1", "RUN5"],
            ),
            (
                pl.concat(
                    [test_master_frame.slice(0, 1), test_master_frame.slice(4, 1)]
                ),
                pl.concat(
                    [test_joint_qc_frame.slice(0, 1), test_joint_qc_frame.slice(4, 1)]
                ),
            ),
        ),
        (
            # an explicit run ID matching nothing selects no rows
            (
                test_master_frame,
                test_joint_qc_frame,
                "dragen",
                None,
                ["RUN_DOES_NOT_EXIST"],
            ),
            (
                test_master_frame.filter(pl.col("WORKFLOW_TYPE") == "NA"),
                test_joint_qc_frame.filter(pl.col("WORKFLOW_TYPE") == "NA"),
            ),
        ),
        (
            # no selector defaults to the last DEFAULT_PLOT_LAST_RUNS (10) runs
            (
                test_master_frame.clone()
                .replace_column(2, pl.Series("WORKFLOW_TYPE", ["dragen"] * 12))
                .replace_column(3, pl.Series("WORKFLOW_VERSION", ["2.6.2.4"] * 12)),
                test_joint_qc_frame.clone()
                .replace_column(1, pl.Series("WORKFLOW_TYPE", ["dragen"] * 12))
                .replace_column(2, pl.Series("WORKFLOW_VERSION", ["2.6.2.4"] * 12)),
                "dragen",
                None,
                None,
            ),
            (
                test_master_frame.slice(2, 10)
                .replace_column(2, pl.Series("WORKFLOW_TYPE", ["dragen"] * 10))
                .replace_column(3, pl.Series("WORKFLOW_VERSION", ["2.6.2.4"] * 10)),
                test_joint_qc_frame.slice(2, 10)
                .replace_column(1, pl.Series("WORKFLOW_TYPE", ["dragen"] * 10))
                .replace_column(2, pl.Series("WORKFLOW_VERSION", ["2.6.2.4"] * 10)),
            ),
        ),
        (
            # Shared run IDs across workflows do not leak between master and joint QC
            (
                test_master_frame.slice(0, 2).replace_column(
                    1, pl.Series("RUN", ["RUN1"] * 2)
                ),
                test_joint_qc_frame.slice(0, 2).replace_column(
                    0, pl.Series("RUN_ID", ["RUN1"] * 2)
                ),
                "dragen",
                None,
                [
                    "RUN1",
                ],
            ),
            (
                test_master_frame.slice(0, 1),
                test_joint_qc_frame.slice(0, 1),
            ),
        ),
    ],
)
def test_select_plot_data(inputs, wants):
    """Master and joint-QC rows are both filtered to the same selected runs."""
    got_0, got_1 = MetricPlots.__new__(MetricPlots).select_plot_data(
        master=inputs[0],
        joint_qc=inputs[1],
        workflow_type=inputs[2],
        plot_last_runs=inputs[3],
        plot_run_ids=inputs[4],
    )
    assert got_0.equals(wants[0])
    assert got_1.equals(wants[1])


@pytest.mark.parametrize(
    ("samples", "samplesheet", "want"),
    [
        (
            # RECORD_TYPE comes from the sample sheet, not metric content or SAMPLE_ID text
            pl.DataFrame(
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
            pl.DataFrame(
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
            pl.DataFrame(
                {
                    "SAMPLE_ID": ["NOT_IN_SAMPLESHEET"],
                }
            ),
            pl.DataFrame(
                {
                    "Sample_ID": ["OTHER_SAMPLE"],
                    "Sample_Type": ["DNA"],
                }
            ),
            ["SAMPLE"],
        ),
        (
            # a sample sheet without a Sample_Type column falls back to unknown
            pl.DataFrame(
                {
                    "SAMPLE_ID": ["SAMPLE01"],
                }
            ),
            pl.DataFrame(
                {
                    "Sample_ID": ["SAMPLE01"],
                }
            ),
            ["SAMPLE"],
        ),
        (
            # a Sample_ID mapped to more than one Sample_Type cannot be classified unambiguously
            pl.DataFrame(
                {
                    "SAMPLE_ID": [
                        "DUPLICATE_ID",
                        "SOLO_SAMPLE",
                    ],
                }
            ),
            pl.DataFrame(
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
            pl.DataFrame(
                {
                    "SAMPLE_ID": ["SAMPLE01"],
                }
            ),
            pl.DataFrame(
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
    got = MetricPlots.__new__(MetricPlots)._add_record_type(samples, samplesheet)
    assert got["RECORD_TYPE"].to_list() == want
