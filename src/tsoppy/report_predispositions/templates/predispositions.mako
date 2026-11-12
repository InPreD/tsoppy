# [${sample_id}] Pipeline: ${tsoppy_version}
# [${sample_id}] Workflow id: ${workflow_id}
# [${sample_id}] Cancer susceptibility genes are defined in:
# [${sample_id}]    file: ${cancer_susceptibility_genes}
# [${sample_id}]    source: ${source}
# [${sample_id}] Small variant calls are defined in: ${small_variant_calls}
# [${sample_id}] Cumulative length of all the targeted coding regions (in millions of bases, seq. depth >= 50): ${ "{:.2f}".format(coding_region_size_in_megabases) }
# [${sample_id}] Tumor purity (as a fraction between 0 and 1): ${ "{:.2f}".format(tumor_purity) }
# [${sample_id}]    The tumor purity is provided as an input parameter to the pipeline.
# [${sample_id}] Column format of gene_predispositions:
# [${sample_id}]    a list of character strings of the form:
# [${sample_id}]    [${ "]_[".join(column_list) }]
# [${sample_id}]    where:
% for column in sorted(all_columns_values):
# [${sample_id}]        ${column}: [ ${ " | ".join(all_columns_values[column]) } ]
% endfor
${ "\t".join(predispositions.columns) }
% for row in predispositions.iter_rows():
${ "\t".join(row) }
% endfor
