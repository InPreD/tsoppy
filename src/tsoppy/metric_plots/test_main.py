from pathlib import Path

from tsoppy.metric_plots.main import MetricPlots


def main():
    metric_plotter = MetricPlots(
        input_directory=Path("tests/test_data/metric_plots/in/"),
        run_id_file=Path("tests/test_data/metric_plots/in/run_ids.txt"),
        output_directory=Path("tests/test_data/metric_plots/out/"),
        create_plots=False,
    )

    metric_plotter.run()

    print("Metric plotting test completed.")
    print("Created:")
    print("out/intermediate_metrics_files/master_metrics_table.tsv")
    print("out/intermediate_metrics_files/joint_sequencing_QC_file.tsv")


if __name__ == "__main__":
    main()
