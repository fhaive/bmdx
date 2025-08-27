# bmdx_package
This repository contain the bmdx package for dose dependent analysis of omics data.
[BMDx paper](https://doi.org/10.1093/bioinformatics/btaa030)

# Manual
The R package manual is available here: [Manual](https://gitlab.com/grecolab_group/bmdx_package/-/blob/master/manual/bmdx_2.0.pdf?ref_type=head)
It was created with the following command: ```devtools::build_manual(path = ".")```

For more informations download the following [html file](https://gitlab.com/grecolab_group/bmdx_package/-/blob/master/manual/manual.html?ref_type=heads)

# Dependencies
```
Depends: R (>= 3.5.0)
Imports: 
    car,
    dendextend,
    dplyr,
    drc,
    ggplot2,
    heatmaply,
    investr,
    limma,
    modelr,
    readxl,
    trend,
    viridis,
    utils,
    minpack.lm,
    stats,
    forcats,
    plotly,
    UpSetR,
    igraph,
    visNetwork,
    networkD3
```

# Installation
Download the source file here: [tar.gz file](https://gitlab.com/grecolab_group/bmdx_package/-/blob/master/bmdx_2.0.tar.gz?ref_type=heads).
The source file was created with the following command: ```devtools::build()```
Use the following command to install the package: ```install.packages("path_to_the_tar.gz.file", type = "source", repos = NULL)```

# Sample Data
Sample data for bleomycin exposure at multiple doses and multiple time points can be found [here](https://gitlab.com/grecolab_group/bmdx_package/-/tree/master/bleomycin_data?ref_type=heads). 

## Authors and acknowledgment
Contributing authors: Michele Fratello, Alisa Pavel

## License
For open source projects, say how it is licensed.

## Project status
This is an ongoing version of the project. If you experiences issues or find any bugs please contact me.

# Example of analysis

```

library(bmdx)
# READ DATA
metadata = read_excel_allsheets(filename =  "bleomycin_data/phenodata_all_timepoints.xlsx")
experimental_data = read_excel_allsheets(filename =  "bleomycin_data/vst_counts_all_timepoints.xlsx", first_col_as_rownames = TRUE,is_rnaseq_raw_count = FALSE)


# CONVERT DATA FOR MODELLING
x = "dose"
y = "expr"
sample_id_col = "sample_id"
dose_id_col = "condition"
other_variables_id_col = c("timepoint")

data_dictionary = create_data_structure(experimental_data,
                                             metadata,
                                             sample_id_col = sample_id_col,
                                             dose_id_col = dose_id_col,
                                             other_variables_id_col = other_variables_id_col ,
                                             x = x,
                                             y = y)

# FILTER BY FOLD CHANGE
samples_variable = "sample_id"
dose_variable = "condition"
time_point_variable = c("timepoint")
pval_th = 0.05
fc_based_filtering = perform_differential_expression_analysis_filtering(data_dictionary,
                                experimental_data,
                                metadata,
                                time_point_variable,
                                dose_variable,
                                samples_variable,
                                fcAdjustment = "fdr",
                                fcPval.th = pval_th,
                                fc.th = 2,
                                nCores = 1,
                                x = "dose",
                                y = "expr",
                                other_variables_id_col = NULL)

plot_filtering_pie_chart(fc_based_filtering$fc_res_unfiltered_dataframe,Pval.th = pval_th)
data_dictionary = fc_based_filtering$filtered_data_dictionary

# BUILD LIST OF MODELS FOR FITTING
model_list = build_models(model_names = c("linear","hill","power",
                                            "poly2","poly3","poly4","poly5",
                                            "exp2","exp3","exp4","exp5",
                                            "llog2","llog3","llog4","llog5",
                                            "mm2","weibul12","weibul13","weibul14",
                                            "weibul22","weibul23","weibul24"), max_iter = 1024, data_type = "continuous", x = x, y = y)

deviation_type = "standard"
rl = 1.349
confidence_interval = 0.95
variance_type = "constant"
significance_level = 0.05
nCores = 1
is_parallel = FALSE

all_fitted_models = fitting_list(data_dictionary,
                                 model_list,
                                 deviation_type = deviation_type,
                                 rl = rl,
                                 confidence_interval = confidence_interval,
                                 variance_type = variance_type,
                                 significance_level = significance_level,
                                 nCores = nCores,
                                 is_parallel = is_parallel)

all_stats = compute_model_statistics(fitted_models = all_fitted_models,
                                      other_variables_id_col = other_variables_id_col,
                                      nCores = 1)

loofth = lower_bound_th = upper_bound_th =  0.1
bmd_bmdl_th = bmdu_bmd_th = 20
bmdu_bmdl_th = 40
r2_th = 0.6
filter_by_monotonicity = F
bmd_na_filter = bmdl_na_filter = bmdu_na_filter = ic50_na_filter = TRUE
filter_by_lack_of_fit = T
r2_filter = T
ratio_filter = filter_lower_bound = filter_upper_bound  = F

filtered_models = model_filtering(fitted_models = all_fitted_models,
                                       loofth = loofth,
                                       lower_bound_th = lower_bound_th,
                                       upper_bound_th = upper_bound_th,
                                       bmd_bmdl_th = bmd_bmdl_th,
                                       bmdu_bmd_th = bmdu_bmd_th,
                                       bmdu_bmdl_th = bmdu_bmdl_th,
                                       filter_lower_bound = filter_lower_bound,
                                       filter_upper_bound = filter_upper_bound,
                                       filter_by_lack_of_fit = filter_by_lack_of_fit,
                                       ratio_filter = ratio_filter,
                                       bmd_na_filter = bmd_na_filter,
                                       bmdl_na_filter = bmdl_na_filter,
                                       bmdu_na_filter = bmdu_na_filter,
                                       ic50_na_filter = ic50_na_filter,
                                       r2_filter = r2_filter,
                                       r2_th = r2_th,
                                       filter_by_monotonicity = filter_by_monotonicity)


all_stats_filtered = compute_model_statistics(filtered_models,
                                              other_variables_id_col = other_variables_id_col,
                                              nCores = 1)

filtered_models_with_avg = add_average_models(filtered_models)

all_stats_filtered_with_avg = compute_model_statistics(filtered_models_with_avg,
                                                       other_variables_id_col = other_variables_id_col,
                                                       nCores = 1)


res = select_optimal_models(filtered_models,
                            method = "AIC",
                            time_col_id = "timepoint",
                            optional_col_ids = NULL,
                            nCores = 1)

optimal_models = res$optimal_models
optimal_models_stats = res$BMD_tab_optimal

#### PLOTTING

library(dplyr)
plot_bmdx(optimal_models[[1]][[1]])

current_stat = optimal_models_stats

p = plot_histogram(current_stat, y_val = "BMD",
                   color_by = "Model",
                   group_by = "timepoint",
                   group_by2 = NULL,
                   filter_column = NULL,
                   filter_by = list(c("linear","poly2","poly3"),c("experiment_name")))
p

current_stat2 = current_stat
current_stat2$timepoint  = as.factor(current_stat2$timepoint)

p=plot_scatter(current_stat2, x_val = "BMDL",y_val = "BMD",
               color_by = "Model",
               group_by = "timepoint",
               group_by2 = NULL,
               filter_column = NULL,
               filter_by = list(c("experiment_name")))
p

p=plot_pie_chart(current_stat, category = "Model",
                 group_by = "timepoint",
                 group_by2 = "Experiment",
                 filter_column = NULL,
                 filter_by = list(c("experiment_name")))
p

p = aggregate_rows_time(current_stat,
                        gen_feat = "Feature",
                        first_feat = "timepoint",
                        group_by = NULL,
                        filter_column = NULL,
                        filter_by = NULL)
p

p = ecdf_plots (mod_stats = current_stat,
                rel_variable = "timepoint",
                group_by = NULL,
                is_group_by_numeric = FALSE,
                other_variables = NULL,
                number_of_column = 2,
                scaling = TRUE, # to use when compounds with different doses are analysed
                filter_column = NULL,
                filter_by = list(c("linear")),
                plot_type = "ecdf")
p

p =  upset_plot(mod_stats = current_stat,
                rel_variable = "timepoint",
                group_by = NULL,
                other_variables = NULL,
                filter_column = NULL,
                filter_by = list(c("linear")),
                nintersects = 3,
                group.by = "degree",
                order.by ="degree")
p

result = compute_gene_frequencies(mod_stats = current_stat,
                                  th = 0.7,
                                  rel_variable = "timepoint",
                                  group_by = "None",
                                  split_by = "None")

result$lollipol_plot_list
gsea_most_freq_genes = gprofiler2::gost(query = result$gene_list,organism = "hsapiens",multi_query = FALSE,ordered_query = TRUE, correction_method = "fdr", domain_scope = "annotated")
gprofiler2::gostplot(gsea_most_freq_genes)

gene_pair_res = gene_pairs_analysis(select_experiment = "bleomycin",
                                     select_time = "24",
                                     filtered_optimal_models = optimal_models,
                                     BMD_tab = optimal_models_stats,
                                     length_vectors = 1000,
                                     nCores = 2,
                                     phenoList = metadata,
                                     doseColID = "condition",
                                     timeColID = "timepoint",
                                     other_variables_id_col = NULL)

stat_table = gene_pair_res$list_gene_pairs_statistics
newdata = gene_pair_res$newdata
selectedrowindex = 1

stringa = rownames(stat_table)[selectedrowindex]
models = unlist(strsplit(stringa, split = "____"))

p = plot_gene_pairs(optimal_models[[models[1]]][[stat_table[selectedrowindex,"Model 1"]]],
                    optimal_models[[models[2]]][[stat_table[selectedrowindex,"Model 2"]]],
                    newdata,
                    main = stat_table[selectedrowindex,1],
                    # main = paste(stat_table[selectedrowindex,"Feature 1"],
                    #              stat_table[selectedrowindex,"Feature 2"],
                    #              sep = " vs "),
                    feat1 = stat_table[selectedrowindex,"Feature 1"],
                    feat2 = stat_table[selectedrowindex,"Feature 2"]

)
p
View(gene_pair_res$list_gene_pairs_statistics)
library(dendextend)
res_gene_clusters = cluster_genes_pairs(gene_pair_res$list_gene_pairs_statistics,2)
res_gene_clusters$heatmap_analysis
res_gene_clusters$clusters

data("ens_human_gene_graph")
library(igraph)
res = make_d3_network(stat_table,g=ens_human_gene_graph, positive = TRUE,th = 0.99)
```

