# BMDx R Package

This repository contains the code for the **bmdx** R package, designed for **dose-dependent analysis of omics data**.  
It is a reimplementation of the code originally provided in the **bmdx Shiny app**.  

📄 [Original BMDx paper](https://doi.org/10.1093/bioinformatics/btaa030)


## 📖 Manual
The package manual is available here:  
[**manual**](https://github.com/fhaive/bmdx/blob/main/manual.0.pdf)

The manual was generated using:
```r
devtools::build_manual(path = ".")
```

---

## 💻 Installation

1. **Create tar.gz file**:  

   ```r
   devtools::build()
   ```

2. **Install the package**:
   ```r
   install.packages("path_to_the_tar.gz.file", type = "source", repos = NULL)
   ```

---

## 📦 Dependencies

The package requires **R (>= 3.5.0)** and the following dependencies:

```
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

---

## 👩‍🔬 Authors & Acknowledgments
- **Angela Serra**  
- **Michele Fratello**  
- **Giorgia Migliaccio**

---

## 📌 Project Status
The project is **ongoing**.

---

## 📊 Sample Data
Example data for bleomycin exposure at multiple doses and time points is available here:  
[**Bleomycin Data**](https://github.com/fhaive/bmdx/tree/main/bleomycin_data)

## 📝 Example of BMDx package Usage

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
other_variables_id_col = NULL

data_dictionary = create_data_structure(experimental_data,
                                             metadata,
                                             sample_id_col = sample_id_col,
                                             dose_id_col = dose_id_col,
                                             other_variables_id_col = other_variables_id_col ,
                                             x = x,
                                             y = y)

# CONVERT DATA FOR MODELLING
x = "dose"
y = "expr"
sample_id_col = "sample_id"
dose_id_col = "condition"
time_col_id = "timepoint"
other_variables_id_col = NULL

data_dictionary = create_data_structure(experimental_data,
                                        metadata,
                                        sample_id_col = sample_id_col,
                                        dose_id_col = dose_id_col,
                                        other_variables_id_col = c(time_col_id,other_variables_id_col) ,
                                        x = x,
                                        y = y)

# FILTER BY FOLD CHANGE
samples_variable = "sample_id"
dose_variable = "condition"
time_point_variable = c("timepoint")
pval_th = 0.05
fc.th =1.5

fc_based_filtering = perform_differential_expression_analysis_filtering(data_dictionary,
                                                                        experimental_data,
                                                                        metadata,
                                                                        time_point_variable,
                                                                        dose_variable,
                                                                        samples_variable,
                                                                        fcAdjustment = "fdr",
                                                                        fcPval.th = pval_th,
                                                                        fc.th = fc.th,
                                                                        nCores = 1,
                                                                        x = "dose",
                                                                        y = "expr",
                                                                        other_variables_id_col = other_variables_id_col)

plot_filtering_pie_chart(tab_unfiltered = fc_based_filtering$fc_res_unfiltered_dataframe,Pval.th = pval_th, title = "Differential analysis")
data_dictionary = fc_based_filtering$filtered_data_dictionary


# BUILD LIST OF MODELS FOR FITTING
model_list = build_models(model_names = c("linear","poly2","hill","exp2","power"), max_iter = 1024, data_type = "continuous", x = x, y = y)

deviation_type = "standard"
rl = 1.349
variance_type = "constant"
confidence_interval = 0.95
significance_level = 0.05 # change the description in the manual
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
                                     other_variables_id_col = c(time_col_id,other_variables_id_col),
                                     nCores = 1,is_parallel = F)

loofth = lower_bound_th = upper_bound_th =  0.1
bmd_bmdl_th = bmdu_bmd_th = 20
bmdu_bmdl_th = 40
r2_th = 0.6
filter_by_monotonicity = F
bmd_na_filter = bmdl_na_filter = bmdu_na_filter = T
ic50_na_filter = F
filter_by_lack_of_fit = F
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
                                              other_variables_id_col = c(time_col_id,other_variables_id_col),
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
current_stat$timepoint = as.factor(current_stat$timepoint)
p = plot_histogram(current_stat, y_val = "BMD",
                   color_by = "timepoint",
                   group_by = "timepoint",
                   group_by2 = NULL,
                   filter_column = NULL,
                   filter_by = list(c("24")))+ ggplot2::theme(
                     axis.title.x = ggplot2::element_text(size = 16),
                     axis.title.y = ggplot2::element_text(size = 16),
                     axis.text.x = ggplot2::element_text(size = 10),
                     axis.text.y = ggplot2::element_text(size = 10),
                     strip.text = ggplot2::element_text(size = 16),  # group_by label size
                     legend.title = ggplot2::element_text(size = 16), # legend title
                     legend.text = ggplot2::element_text(size = 14) 
p

p=plot_pie_chart(current_stat, category = "Model",
                 group_by = "timepoint",
                 group_by2 = NULL,
                 filter_column = NULL,
                 filter_by = list(c("24")))+ ggplot2::theme(
                   axis.title.x = ggplot2::element_text(size = 16),
                   axis.title.y = ggplot2::element_text(size = 16),
                   axis.text.x = ggplot2::element_text(size = 10),
                   axis.text.y = ggplot2::element_text(size = 10),
                   strip.text = ggplot2::element_text(size = 16),  # group_by label size
                   legend.title = ggplot2::element_text(size = 16), # legend title
                   legend.text = ggplot2::element_text(size = 14) 
                 )
p

p = aggregate_rows_time(current_stat,
                        gen_feat = "Feature",
                        first_feat = "timepoint",
                        group_by = NULL,
                        filter_column = NULL,
                        filter_by = NULL) + ggplot2::theme(
                          axis.title.x = ggplot2::element_text(size = 16),
                          axis.title.y = ggplot2::element_text(size = 16),
                          axis.text.x = ggplot2::element_text(size = 16),
                          axis.text.y = ggplot2::element_text(size = 16)
                        )
p

p=plot_scatter(current_stat, x_val = "BMDL",y_val = "BMD",
               color_by = "Model",
               group_by = "timepoint",
               group_by2 = NULL,
               filter_column = NULL,
               filter_by = list(c("experiment_name")))
p

p = ecdf_plots (mod_stats = current_stat,
                rel_variable = "timepoint",
                group_by = NULL,
                is_group_by_numeric = FALSE,
                other_variables = NULL,
                number_of_column = 2,
                scaling = F, # to use when compounds with different doses are analysed
                filter_column = NULL,
                plot_type = "ecdf",linewidth = 1.2)+ ggplot2::theme(
                  axis.title.x = ggplot2::element_text(size = 16),
                  axis.title.y = ggplot2::element_text(size = 16),
                  axis.text.x = ggplot2::element_text(size = 16),
                  axis.text.y = ggplot2::element_text(size = 16)
                )
p

p =  upset_plot(mod_stats = current_stat,
                rel_variable = "timepoint",
                group_by = NULL,
                other_variables = NULL,
                filter_column = NULL,
                filter_by = list(),
                nintersects = 10,
                group.by = "degree",
                order.by ="degree",text.scale = 2)
p

#TPOD
library(scam)

# Parameters for tPOD computation
pod_value <- "BMD" #or "BMDL", "BMDU"
percentile <- 0.20 # a number between 0 and 1
lowest_method = "lowest" #or "LCRD"

model_stats_list <- setNames(
  split(optimal_models_stats, interaction(optimal_models_stats[["Experiment"]], 
                                          optimal_models_stats[[time_col_id]], 
                                          drop = TRUE, sep = "_")),
  levels(interaction(optimal_models_stats[["Experiment"]], 
                     optimal_models_stats[[time_col_id]], 
                     drop = TRUE, sep = "_"))
)

tpod_methods_list <- c("percentile", "first_mode", "lowest", "accumulation")

# Compute tPOD for each experimental condition and method
result_tPOD <- lapply(model_stats_list, function(model_stats) {
  res = apply_tpod_methods(model_stats = model_stats,
                     tpod_methods_list = tpod_methods_list,
                     pod_value = "BMD",
                     percentile = percentile,
                     lowest_method = lowest_method)
  rownames(res) = res$Method
  return(res)
})


p1 = tpod_plot(pod_vector = model_stats_list$bleomycin_24$BMD,
          tpod_method = "accumulation",
          pod_value = "BMD",
          tPOD = as.numeric(result_tPOD$bleomycin_24["accumulation","tPOD"]),
          xlog = T,
          subtitle = "tPOD bleomycin exposure 24h") + theme(legend.position = "bottom")


BMD_values <- model_stats_list$bleomycin_24$BMD

lowest = as.numeric(result_tPOD$bleomycin_24["lowest","tPOD"])
percentile = as.numeric(result_tPOD$bleomycin_24["percentile","tPOD"])
mean_value = as.numeric(result_tPOD$bleomycin_24["mean","tPOD"])
accumulation = as.numeric(result_tPOD$bleomycin_24["accumulation","tPOD"])
first_mode = as.numeric(result_tPOD$bleomycin_24["first_mode","tPOD"])

p = plot_BMD_tPOD_density(BMD_values, lowest = lowest,
                      percentile = percentile,
                      mean_value = mean_value,
                      accumulation = accumulation,
                      first_mode = first_mode)

# Gene frequency

result = compute_gene_frequencies(mod_stats = current_stat,
                                  th = 0.7,
                                  rel_variable = "timepoint",
                                  group_by = "None",
                                  split_by = "None")

result$lollipol_plot_list
gsea_most_freq_genes = gprofiler2::gost(query = result$gene_list,organism = "hsapiens",multi_query = FALSE,ordered_query = TRUE, correction_method = "fdr", domain_scope = "annotated")
gprofiler2::gostplot(gsea_most_freq_genes)

# Gene pairs

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
res_gene_clusters = cluster_genes_pairs(comparison_pairs = gene_pair_res$list_gene_pairs_statistics,
                                        nclust = 3,method = "correlation")

res_gene_clusters$heatmap_analysis
res_gene_clusters$complex_heatmap

complex_heatmap = res_gene_clusters$complex_heatmap
complex_heatmap@matrix_param$show_row_names <- FALSE
complex_heatmap@matrix_param$show_column_names <- FALSE
draw(complex_heatmap)

res_gene_clusters$n_genes_par_cluster
clusters = res_gene_clusters$clusters

min_n_elem = 20
tab = table(clusters)
to_rem = which(tab<min_n_elem)

if(length(to_rem)>0){
  clusters = clusters[(clusters %in% to_rem)==FALSE]
}


query  =  list()
for(i in min(clusters):max(clusters)){
  query[[paste("Cluster",i,sep = "_")]]  =  names(clusters)[clusters  ==  i]
}

organism  =  organism = "hsapiens"
correction_method  =  correction_method = "fdr"

res_GSEA_over_corr_clust  =  gprofiler2::gost(query  =  query,
                                              organism  =  organism,
                                              ordered_query  =  FALSE,
                                              correction_method  =  correction_method,
                                              domain_scope  =  "annotated")

res_GSEA_over_corr_clust$result$query <- as.factor(res_GSEA_over_corr_clust$result$query)
p <- plot_enrichment_gost_bubble(data = res_GSEA_over_corr_clust$result,
                                 source_filter = "GO:BP",
                                 top_n = 10,
                                 facet_column = "query")

p
```

