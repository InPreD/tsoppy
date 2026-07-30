# in the folder in/both* I have added some files from dragen and localapp to see if code can correctly process them
# the master metrics table and joint sequencing QC file should be generated with both files.

tsoppy metric-plots --input-directory tests/test_data/metric_plots/in/both_dragen_localapp/
--run-id-file tests/test_data/metric_plots/in/run_ids.txt --output-directory tests/test_data/metric_plots/out/ --no-create-plots