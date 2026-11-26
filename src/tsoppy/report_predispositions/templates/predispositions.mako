# [${sample_id}] Pipeline: ${tsoppy_version}
# [${sample_id}] Workflow id: ${workflow_id}
# [${sample_id}] Cancer susceptibility genes are defined in:
# [${sample_id}]    file: ${cancer_susceptibility_genes}
# [${sample_id}]    source: ${source}
% if small_variant_calls_normal_DNA is not None:
# [${sample_id}] Small variant calls of normal DNA sample are defined in: ${small_variant_calls_normal_DNA}
% endif
% if small_variant_calls_tumor_DNA is not None:
# [${sample_id}] Small variant calls of tumor DNA sample are defined in: ${small_variant_calls_tumor_DNA}
% endif
% if sample_id_tumor_RNA is not None:
# [${sample_id}] Tumor RNA sample id: ${sample_id_tumor_RNA}
% endif
# [${sample_id}] Cumulative length of all the targeted coding regions (in millions of bases, seq. depth >= 50): ${ "{:.2f}".format(coding_region_size_in_megabases) }
% if small_variant_calls_tumor_DNA is not None:
# [${sample_id}] Tumor purity (as a fraction between 0 and 1): ${ "{:.2f}".format(tumor_purity) }
# [${sample_id}]    The tumor purity is provided as an input parameter to the pipeline.
% endif
# [${sample_id}] Column format of gene_predispositions:
# [${sample_id}]    a list of character strings of the form:
# [${sample_id}]    [${ "]_[".join(column_list) }]
# [${sample_id}]    where:
% for column in sorted(all_columns_values):
# [${sample_id}]        ${column}: [ ${ " | ".join(all_columns_values[column]) } ]
% endfor
${ "\t".join(predispositions.columns) }
% if predispositions.is_empty():
# [${sample_id}] No variants found.
% else:
    % for row in predispositions.iter_rows():
${ "\t".join(row) }
    % endfor
