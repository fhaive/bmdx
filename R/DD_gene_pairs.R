#' Perform gene pairs analysis for selected experiments and times
#' This function performs gene pairs analysis for selected experiments and times,
#' using the provided filtered_optimal_models, BMD_tab, length_vectors, nCores,
#' phenoList, doseColID, timeColID, and other_variables_id_col.
#' @param select_experiment A character vector specifying the selected
#' experiments for the analysis.
#' @param select_time A numeric vector specifying the selected times for the analysis.
#' @param filtered_optimal_models A list of filtered optimal models.
#' @param BMD_tab The BMD (Benchmark Dose) table.
#' @param length_vectors An integer specifying the length of vectors.
#' @param nCores An integer specifying the number of cores to use for parallel processing.
#' @param phenoList A list containing the phenotype data.
#' @param doseColID The column ID for the dose in the phenotype data.
#' @param timeColID The column ID for the time in the phenotype data.
#' @param other_variables_id_col The column ID for other variables in the phenotype data.
#' @return A list containing gene pairs statistics and the newdata used for the analysis.
#' @export

gene_pairs_analysis = function(select_experiment,
                               select_time,
                               filtered_optimal_models,
                               BMD_tab,
                               length_vectors,
                               nCores,
                               phenoList,
                               doseColID,
                               timeColID,
                               other_variables_id_col){

  phenoDt <- phenoList[[1]]

  newdata = data.frame(dose = seq(min(phenoDt[,doseColID]), max(phenoDt[,doseColID]), length.out = length_vectors))

  experiment = unlist(lapply(strsplit(x = names(filtered_optimal_models),split = "_"),function(elem)elem[1]))
  times = unlist(lapply(strsplit(x = names(filtered_optimal_models),split = "_"),function(elem)elem[2]))

  # combo = unique(expand.grid(experiment, times))
  combo = unique(expand.grid(select_experiment,select_time))
  list_gene_pairs_statistics = c()

  for (ii in 1:nrow(combo)) {
    idx = which(times == combo[ii,2] & experiment == combo[ii,1])

    if (length(idx) > 0) {
      models_tp3 = filtered_optimal_models[idx]
      for (i in 1:length(models_tp3)) {
        models_tp3[[i]] = find_best_model_aic(models_tp3[[i]])
      }
      statistics_gene_pairs = gene_pairs_comparison(filtered_optimal_models = models_tp3,
                                                    other_variables_id_col = other_variables_id_col,
                                                    newdata,
                                                    nCores = 1)
      # list_gene_pairs_statistics[[paste(combo[ii,1],combo[ii,2],sep="_")]] = statistics_gene_pairs
      list_gene_pairs_statistics = rbind(list_gene_pairs_statistics, statistics_gene_pairs)
    }

  }

  results = list("list_gene_pairs_statistics" = list_gene_pairs_statistics,
                 "newdata" = newdata)
  return(results)
}

message_parallel <- function(...){
  system(sprintf('echo "\n%s\n"', paste0(..., collapse = "")))
}

#' Compare gene pairs for their dose-dependent patterns
#' This function compares gene pairs for their dose-dependent patterns using the
#' provided filtered_optimal_models, other_variables_id_col, and newdata.
#' It calculates the difference, behavior, correlation, and coefficient for each
#' pair of genes in the models.
#' @param filtered_optimal_models A list of filtered optimal models containing
#' gene-related information.
#' @param other_variables_id_col The column ID for other variables in the gene data.
#' @param newdata A data frame containing dose information for new observations.
#' @param nCores An integer specifying the number of cores to use for parallel
#' processing (default is 40).
#' @return A data frame containing gene pairs statistics, including Experiment 1,
#' Model 1, Experiment 2, Model 2, Difference Trend, Coefficient, and Correlation of Gene Patterns.
#' @export

gene_pairs_comparison = function(filtered_optimal_models,
                                 other_variables_id_col,
                                 newdata,nCores = 40){

  res = parallel::mclapply(X = 1:(length(filtered_optimal_models)), FUN = function(i){
    message_parallel(paste(i,"/",length(filtered_optimal_models),sep = ""))

    opt_model_i = filtered_optimal_models[[i]]
    elem = unlist(strsplit(names(filtered_optimal_models)[i],split = "_"))
    elem_i = elem#[length(elem)]
    yhat_i = predict(opt_model_i$fitted, newdata = newdata)

    names_set = rep("",length((i + 1):length(filtered_optimal_models)))
    df_statistics = c()#matrix(0, nrow = length((i+1):length(filtered_optimal_models)), ncol = 7)

    index = 1

    for (j in 1:length(filtered_optimal_models)) {
      opt_model_j = filtered_optimal_models[[j]]
      elem = unlist(strsplit(names(filtered_optimal_models)[j],split = "_"))
      elem_j = elem#[length(elem)]

      yhat_j = predict(opt_model_j$fitted, newdata = newdata)
      diff = abs(yhat_i - yhat_j)
      df = data.frame(x = newdata$dose,y = diff)

      fitlm = stats::lm(y~x,data = df)
      behaviour = sign(fitlm$coefficients[2]) # is the difference increasing or decreasing with time

      pcor = stats::cor(yhat_i ,yhat_j) # correlation between the two dose-dependent gene profile patterns

      normalized_euclidean <- function(a, b) {
        d <- sqrt(sum((a - b)^2))
        d / sqrt(length(a))  # Normalized by sqrt of vector length
      }

      n_eucl_dist = normalized_euclidean(yhat_i ,yhat_j)
      # df_statistics[index,] = c(elem_i, opt_model_i$mod_name, elem_j, opt_model_j$mod_name, behaviour, fitlm$coefficients[2], pcor)
      df_statistics = rbind(df_statistics,c(elem_i, opt_model_i$mod_name,
                                            elem_j, opt_model_j$mod_name,
                                            behaviour, fitlm$coefficients[2], pcor, n_eucl_dist))

      names_set[index] =  paste(names(filtered_optimal_models)[i],names(filtered_optimal_models)[j],sep = "____")

      index = index + 1

    }

    colnames(df_statistics) = c("Experiment 1",paste(other_variables_id_col, "1", sep = " "),"Feature 1","Model 1",
                                "Experiment 2",paste(other_variables_id_col, "2", sep = " "),"Feature 2","Model 2",
                                "DiffTrend","Coefficient","CorGenePatteerns", "NormalizedEuclideanDistance")


    rownames(df_statistics) = names_set

    return(list(df_statistics))

  },mc.cores = nCores)

  StatMat = do.call(rbind, lapply(X = res, FUN = function(elem)elem[[1]]))
  StatMat = as.data.frame(StatMat)
  StatMat$DiffTrend = as.numeric(as.vector(StatMat$DiffTrend))
  StatMat$Coefficient = as.numeric(as.vector(StatMat$Coefficient))
  StatMat$CorGenePatteerns = as.numeric(as.vector(StatMat$CorGenePatteerns))
  StatMat$NormalizedEuclideanDistance = as.numeric(as.vector(StatMat$NormalizedEuclideanDistance))

  return(StatMat)
}

# gene_pairs_comparison = function(filtered_optimal_models,
#                                  other_variables_id_col,
#                                  newdata,nCores = 40){
#
#   res = parallel::mclapply(X = 1:(length(filtered_optimal_models) - 1), FUN = function(i){
#     message_parallel(paste(i,"/",length(filtered_optimal_models),sep = ""))
#
#     opt_model_i = filtered_optimal_models[[i]]
#     elem = unlist(strsplit(names(filtered_optimal_models)[i],split = "_"))
#     elem_i = elem#[length(elem)]
#     yhat_i = predict(opt_model_i$fitted, newdata = newdata)
#
#     names_set = rep("",length((i + 1):length(filtered_optimal_models)))
#     df_statistics = c()#matrix(0, nrow = length((i+1):length(filtered_optimal_models)), ncol = 7)
#
#     index = 1
#
#     for (j in (i + 1):length(filtered_optimal_models)) {
#       opt_model_j = filtered_optimal_models[[j]]
#       elem = unlist(strsplit(names(filtered_optimal_models)[j],split = "_"))
#       elem_j = elem#[length(elem)]
#
#       yhat_j = predict(opt_model_j$fitted, newdata = newdata)
#       diff = abs(yhat_i - yhat_j)
#       df = data.frame(x = newdata$dose,y = diff)
#
#       fitlm = stats::lm(y~x,data = df)
#       behaviour = sign(fitlm$coefficients[2]) # is the difference increasing or decreasing with time
#
#       pcor = stats::cor(yhat_i ,yhat_j) # correlation between the two dose-dependent gene profile patterns
#
#       # df_statistics[index,] = c(elem_i, opt_model_i$mod_name, elem_j, opt_model_j$mod_name, behaviour, fitlm$coefficients[2], pcor)
#       df_statistics = rbind(df_statistics,c(elem_i, opt_model_i$mod_name,
#                                             elem_j, opt_model_j$mod_name,
#                                             behaviour, fitlm$coefficients[2], pcor))
#
#       names_set[index] =  paste(names(filtered_optimal_models)[i],names(filtered_optimal_models)[j],sep = "____")
#
#       index = index + 1
#
#     }
#
#     colnames(df_statistics) = c("Experiment 1",paste(other_variables_id_col, "1", sep = " "),"Feature 1","Model 1",
#                                 "Experiment 2",paste(other_variables_id_col, "2", sep = " "),"Feature 2","Model 2",
#                                 "DiffTrend","Coefficient","CorGenePatteerns")
#
#
#     rownames(df_statistics) = names_set
#
#     return(list(df_statistics))
#
#   },mc.cores = nCores)
#
#   StatMat = do.call(rbind, lapply(X = res, FUN = function(elem)elem[[1]]))
#   StatMat = as.data.frame(StatMat)
#   StatMat$DiffTrend = as.numeric(as.vector(StatMat$DiffTrend))
#   StatMat$Coefficient = as.numeric(as.vector(StatMat$Coefficient))
#   StatMat$CorGenePatteerns = as.numeric(as.vector(StatMat$CorGenePatteerns))
#
#
#   return(StatMat)
# }

#' Plot gene pairs comparison
#'
#' This function plots the comparison of two gene pairs given the models mod_i
#' and mod_j and newdata containing dose information for new observations.
#' @param mod_i The model for gene pair 1.
#' @param mod_j The model for gene pair 2.
#' @param newdata A data frame containing dose information for new observations.
#' @param main A character string specifying the main title of the plot
#'  (default is an empty string).
#' @param feat1 A character string specifying the label for gene pair 1 in the
#' plot (default is "g1").
#' @param feat2 A character string specifying the label for gene pair 2 in the plot
#' (default is "g2").
#' @return A ggplot object showing the comparison of gene pairs.
#' @export
plot_gene_pairs = function(mod_i,mod_j, newdata, main = "", feat1 = "g1", feat2 = "g2"){#, plot_intervals = FALSE,plot_diff = FALSE){

  yhat_i = predict(mod_i, newdata = newdata)#,interval = "prediction")
  yhat_j = predict(mod_j, newdata = newdata)#,interval = "prediction")
  colnames(yhat_i) = colnames(yhat_j)  = "fit"
  diff = abs(yhat_i - yhat_j)

  df = data.frame(dose = newdata$dose, expression = c(yhat_i[,"fit"],yhat_j[,"fit"]), feature = c(rep(feat1,length(yhat_i[,"fit"])),rep(feat2,length(yhat_j[,"fit"]))))

  p = ggplot2::ggplot(df, ggplot2::aes(x = dose, y = expression, group = feature, color = feature)) +
    ggplot2::geom_line(size = 2) +
    viridis::scale_color_viridis(discrete = TRUE) +
    ggplot2::ggtitle(main) + ggplot2::theme(legend.position = "bottom")

  return(p)
}

#' Load Protein-Protein Interaction (PPI) Graph for an Organism
#'
#' This function loads a precompiled protein-protein interaction (PPI) network graph
#' for a specified organism and gene identifier type (e.g., ENSEMBL, ENTREZ, or SYMBOL).
#' It returns the corresponding `igraph` object representing the network.
#'
#' @param organism A string indicating the organism name. Supported values are
#'        `"human"`, `"mouse"`, and `"rat"`.
#' @param gene_id_type A string indicating the gene identifier type. Supported values are
#'        `"ENSEMBL"`, `"ENTREZ"`, and `"SYMBOL"`.
#' @param graph_type A string indicating the graph type. Supported values are
#'        `"ppi"`, and `"tf"`.
#'
#' @return An `igraph` object containing the PPI network for the specified organism and gene ID type.
#'
#' @details
#' This function expects that the relevant data objects (e.g., `ens_human_gene_graph`,
#' `sym_mouse_gene_graph`, etc.) are available in the environment or included as package data.
#' The function uses `data()` to load the correct graph into the environment before returning it.
#'
#' @export
load_ppi_data = function(organism, gene_id_type, graph_type = "ppi") {


  if(graph_type %in% c("ppi","tf","mirna") == FALSE){
    print("graph type must be one between ppi, tf, and mirna")
    return()
  }

  if(gene_id_type %in% c("ENSEMBL","ENTREZ","SYMBOL") == FALSE){
    print("GeneID type must be one between ENSEMBL,ENTREZ, and SYMBOL")
    return()
  }

  if(organism %in% c("human","mouse","rat") == FALSE){
    print("Organism must be one between human,mouse, and rat")
    return()
  }

  human_ppi = list("ENSEMBL" = "ens_human_gene_graph",
                   "ENTREZID"  = "ent_human_gene_graph",
                   "SYMBOL"  = "sym_human_gene_graph")

  mouse_ppi = list("ENSEMBL" = "ens_mouse_gene_graph",
                   "ENTREZID"  = "ent_mouse_gene_graph",
                   "SYMBOL"  = "sym_mouse_gene_graph")

  rat_ppi   = list("ENSEMBL" = "ens_rat_gene_graph",
                   "ENTREZID"  = "ent_rat_gene_graph",
                   "SYMBOL"  = "sym_rat_gene_graph")

  human_tf = list("ENSEMBL" = "ens_human_gene_graph_tf",
                   "ENTREZID"  = "ent_human_gene_graph_tf",
                   "SYMBOL"  = "sym_human_gene_graph_tf")

  mouse_tf = list("ENSEMBL" = "ens_mouse_gene_graph_tf",
                   "ENTREZID"  = "ent_mouse_gene_graph_tf",
                   "SYMBOL"  = "sym_mouse_gene_graph_tf")

  rat_tf   = list("ENSEMBL" = "ens_rat_gene_graph_tf",
                   "ENTREZID"  = "ent_rat_gene_graph_tf",
                   "SYMBOL"  = "sym_rat_gene_graph_tf")

  human_mirna = list("ENSEMBL" = "ens_human_gene_graph_mirna",
                  "ENTREZID"  = "ent_human_gene_graph_mirna",
                  "SYMBOL"  = "sym_human_gene_graph_mirna")

  mouse_mirna = list("ENSEMBL" = "ens_mouse_gene_graph_mirna",
                  "ENTREZID"  = "ent_mouse_gene_graph_mirna",
                  "SYMBOL"  = "sym_mouse_gene_graph_mirna")

  rat_mirna   = list("ENSEMBL" = "ens_rat_gene_graph_mirna",
                  "ENTREZID"  = "ent_rat_gene_graph_mirna",
                  "SYMBOL"  = "sym_rat_gene_graph_mirna")

  available_graphs = list("ppi" = list(),
                          "tf"  = list(),
                          "mirna" = list())

  available_graphs$ppi[["human"]] = human_ppi
  available_graphs$ppi[["mouse"]] = mouse_ppi
  available_graphs$ppi[["rat"]]   = rat_ppi

  available_graphs$tf[["human"]] = human_tf
  available_graphs$tf[["mouse"]] = mouse_tf
  available_graphs$tf[["rat"]]   = rat_tf

  available_graphs$mirna[["human"]] = human_mirna
  available_graphs$mirna[["mouse"]] = mouse_mirna
  available_graphs$mirna[["rat"]]   = rat_mirna

  dataset_name = available_graphs[[graph_type]][[organism]][[gene_id_type]]
  data(list = dataset_name)

  graph = get(dataset_name)
  # if (organism == "human") {
  #   if (gene_id_type == "ENSEMBL") {
  #     data("ens_human_gene_graph")
  #     graph = ens_human_gene_graph
  #   }
  #   if (gene_id_type == "ENTREZ") {
  #     data("ent_human_gene_graph")
  #     graph = ent_human_gene_graph
  #   }
  #   if (gene_id_type == "SYMBOL") {
  #     data("sym_human_gene_graph")
  #     graph = sym_human_gene_graph
  #   }
  # }
  #
  # if (organism == "mouse") {
  #   if (gene_id_type == "ENSEMBL") {
  #     data("ens_mouse_gene_graph")
  #     graph = ens_mouse_gene_graph
  #   }
  #   if (gene_id_type == "ENTREZ") {
  #     data("ent_mouse_gene_graph")
  #     graph = ent_mouse_gene_graph
  #   }
  #   if (gene_id_type == "SYMBOL") {
  #     data("sym_mouse_gene_graph")
  #     graph = sym_mouse_gene_graph
  #   }
  # }
  #
  # if (organism == "rat") {
  #   if (gene_id_type == "ENSEMBL") {
  #     data("ens_rat_gene_graph")
  #     graph = ens_rat_gene_graph
  #   }
  #   if (gene_id_type == "ENTREZ") {
  #     data("ent_rat_gene_graph")
  #     graph = ent_rat_gene_graph
  #   }
  #   if (gene_id_type == "SYMBOL") {
  #     data("sym_rat_gene_graph")
  #     graph = sym_rat_gene_graph
  #   }
  # }

  return(graph)
}


#' Find the best model based on AIC (Akaike Information Criterion)
#' This function finds the best model based on the Akaike Information Criterion (AIC) among a list of models.
#' @param model_list A list containing multiple models to be compared.
#' @return The model from the list that has the lowest AIC value.
#' @export
find_best_model_aic = function(model_list){
  opt_idx = which.min(unlist(lapply(model_list, function(elem) elem$aic)))
  return(model_list[[opt_idx]])
}

#' Cluster Gene Pairs Based on Expression Profile Similarity
#'
#' This function clusters gene pairs using similarity metrics derived from
#' pairwise comparisons (e.g., correlation and distance). It constructs a matrix
#' representing gene-gene similarity, performs hierarchical clustering, and generates
#' a `ComplexHeatmap` annotated with cluster and centrality information. Optionally, it
#' overlays protein-protein interaction (PPI) data to highlight known gene interactions.
#'
#' @param comparison_pairs A data frame of pairwise gene comparisons. Must contain
#'        columns `"Feature 1"`, `"Feature 2"`, `"CorGenePatteerns"` (correlation), and
#'        `"NormalizedEuclideanDistance"` (distance).
#' @param nclust Integer. Number of clusters to extract from hierarchical clustering. Default is 2.
#' @param method A string specifying the clustering metric: `"euclidean"`, `"correlation"`, or `"combination"`.
#'        If `"combination"`, both normalized distance and inverse correlation are averaged.
#' @param plot_ppi_info Logical. If `TRUE`, adds PPI network information and node centrality measures
#'        (degree, closeness, eigenvector) to the heatmap annotations. Default is `TRUE`.
#' @param gene_id_type A string indicating the type of gene identifier used (e.g., `"ENSEMBL"`, `"ENTREZ"`, `"SYMBOL"`).
#' @param organism A string specifying the organism: `"human"`, `"mouse"`, or `"rat"`.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{clusters}{A named vector of cluster assignments for each gene.}
#'   \item{n_genes_par_cluster}{A table of gene counts per cluster.}
#'   \item{average_correlation_in_clusters}{A numeric vector of mean intra-cluster similarity values.}
#'   \item{list_correlation_in_clusters}{A list of submatrices showing correlation values within each cluster.}
#'   \item{complex_heatmap}{A `ComplexHeatmap` object visualizing the clustered similarity matrix.}
#' }
#'
#' @details
#' The function builds symmetric matrices of pairwise correlation and Euclidean distance,
#' performs hierarchical clustering on the combined or selected metric, and annotates the heatmap
#' with cluster statistics and optional PPI features. When `plot_ppi_info = TRUE`, it highlights
#' interactions using dot markers and colors nodes based on centrality measures derived from organism-specific
#' PPI graphs.
#' @export
cluster_genes_pairs = function(comparison_pairs,
                               nclust=2, method = "euclidean",
                               plot_ppi_info = TRUE,
                               gene_id_type = "ENSEMBL",
                               organism = "human"){



  M = comparison_pairs[,c("Feature 1","Feature 2","CorGenePatteerns")]
  colnames(M) = c("feat1","feat2","cor")
  unique_feat = unique(c(M$feat1, M$feat2))
  M_cor = matrix(0, nrow = length(unique_feat), ncol = length(unique_feat),dimnames = list(unique_feat, unique_feat))
  for (i in 1:nrow(M)) {
    M_cor[M[i,1],M[i,2]] = M[i,3]
    M_cor[M[i,2],M[i,1]] = M[i,3]
  }

  M = comparison_pairs[,c("Feature 1","Feature 2","NormalizedEuclideanDistance")]
  colnames(M) = c("feat1","feat2","NEucDist")
  unique_feat = unique(c(M$feat1, M$feat2))
  M_dist = matrix(0, nrow = length(unique_feat), ncol = length(unique_feat),dimnames = list(unique_feat, unique_feat))
  for (i in 1:nrow(M)) {
    M_dist[M[i,1],M[i,2]] = M[i,3]
    M_dist[M[i,2],M[i,1]] = M[i,3]
  }

  if(method == "euclidean"){
    M2 = M_dist
  }
  if(method == "correlation"){
    M2 = M_cor
  }
  if(method == "combination"){

    M_cor = 1 - M_cor
    M_cor = M_cor/max(M_cor)
    M_dist = M_dist/max(M_dist)
    M2 = (M_cor + M_dist)/2
  }

  row_dend  <- M2 %>%
    dist %>%
    hclust %>%
    as.dendrogram %>%
    set("branches_k_color", k = as.numeric(nclust)) %>%
    ladderize
  # rotate_DendSer(ser_weight = dist(x))
  col_dend  <- M2 %>%
    t %>%
    dist %>%
    hclust %>%
    as.dendrogram %>%
    set("branches_k_color", k = as.numeric(nclust)) %>%
    ladderize

  # heatmap_analysis =  heatmaply::heatmaply(
  #   as.data.frame(M2),
  #   Rowv = row_dend,
  #   Colv = col_dend
  # )

  clusters = dendextend::cutree(col_dend,k = as.numeric(nclust))
  n_genes_par_cluster = table(clusters)

  average_correlation_in_clusters = c()
  for(i in 1:length(unique(clusters))){
    gi = names(clusters)[clusters==i]
    avg_cor = mean(M2[gi,gi])
    average_correlation_in_clusters = c(average_correlation_in_clusters,avg_cor)

  }

  list_correlation_in_clusters = list()
  for(i in 1:length(unique(clusters))){
    gi = names(clusters)[clusters==i]
    list_correlation_in_clusters[[i]]  = M2[gi,gi]
  }


  library(ComplexHeatmap)
  library(circlize)  # For color scales

  # Perform hierarchical clustering
  row_hclust <- hclust(dist(M2))
  col_hclust <- hclust(dist(t(M2)))

  # Convert dendrograms for ComplexHeatmap
  row_dend <- as.dendrogram(row_hclust)
  col_dend <- as.dendrogram(col_hclust)

  # Assign clusters (adjust k as needed)
  row_clusters <- cutree(row_hclust, k = nclust)
  col_clusters <- cutree(col_hclust, k = nclust)

  # Compute cluster sizes
  row_cluster_sizes <- table(row_clusters)
  col_cluster_sizes <- table(col_clusters)

  # Compute mean M2 values per cluster
  compute_cluster_means <- function(matrix_data, cluster_labels) {
    cluster_ids <- unique(cluster_labels)
    cluster_means <- numeric(length(cluster_ids))

    for (i in seq_along(cluster_ids)) {
      cluster <- cluster_ids[i]
      members <- which(cluster_labels == cluster)
      cluster_means[i] <- mean(matrix_data[members, ], na.rm = TRUE)  # Compute mean
    }

    names(cluster_means) <- cluster_ids
    return(cluster_means)
  }

  # Compute row-wise and column-wise mean values
  row_avg_value <- compute_cluster_means(M2, row_clusters)
  col_avg_value <- compute_cluster_means(t(M2), col_clusters)

  # Create row annotations (Cluster Size + Avg Value)
  # row_anno <- rowAnnotation(
  #   Cluster = factor(row_clusters),
  #   Cluster_Size = row_cluster_sizes[row_clusters],
  #   Avg_Value = round(row_avg_value[row_clusters], 2),
  #   col = list(Cluster = structure(1:nclust, names = unique(row_clusters)))
  # )

  # Create column annotations (Cluster Size + Avg Value)
  # col_anno <- HeatmapAnnotation(
  #   Cluster = factor(col_clusters),
  #   Cluster_Size = col_cluster_sizes[col_clusters],
  #   Avg_Value = round(col_avg_value[col_clusters], 2),
  #   col = list(Cluster = structure(1:nclust, names = unique(col_clusters)))
  # )

  col_anno <- HeatmapAnnotation(
    Cluster = factor(col_clusters),
    # Cluster_Size = col_cluster_sizes[col_clusters],
    # Avg_Value = round(col_avg_value[col_clusters], 2),
    col = list(Cluster = structure(1:nclust, names = unique(col_clusters)))
  )

  if(plot_ppi_info){

    graph = load_ppi_data(organism,gene_id_type,graph_type = "ppi")

    edge_df = comparison_pairs[,c("Feature 1","Feature 2")]
    edge_presence = check_edges_in_graph(graph, edge_df)
    colnames(edge_presence) = c("gene1","gene2","presence")

    library(reshape2)
    annotate_matrix <- acast(edge_presence, gene1 ~ gene2,
                             value.var = "presence")
    diag(annotate_matrix) = FALSE
    annotate_matrix[is.na(annotate_matrix)] = FALSE

    graph_tf = load_ppi_data(organism,gene_id_type,graph_type = "tf")
    graph_mirna = load_ppi_data(organism,gene_id_type,graph_type = "mirna")

    # Get heads of all edges
    tf_list    <- unique(ends(graph_tf, E(graph_tf), names = T)[, 1])
    mirna_list <- unique(ends(graph_mirna, E(graph_mirna), names = FALSE)[, 1])


    col_fun_tf <- c("0" = "#9381ff", "1" = "#301934")         # light gray for FALSE, dark violet for TRUE
    col_fun_mirna <- c("0" = "#a8dadc", "1" = "#bc6c25")


    right_anno = rowAnnotation("TF" = as.numeric(rownames(M2) %in% tf_list),
                               "Mirna" = as.numeric(rownames(M2) %in% mirna_list),
                               col = list("TF" = col_fun_tf,
                                          "Mirna" = col_fun_mirna ))

    # Create row annotations (Cluster Size + Avg Value)

    degree_vals = rep(0, nrow(M2))
    names(degree_vals) = rownames(M2)

    closeness_vals = degree_vals
    eigenvector_vals = degree_vals

    vertices = rownames(M2)[which(rownames(M2) %in% V(graph)$name)]

    degree_vals[vertices] = igraph::degree(graph,v = vertices)
    closeness_vals[vertices] = igraph::closeness(graph,v = vertices, normalized = T)
    eigenvector_vals[vertices] = eigen_centrality(graph)$vector[vertices]

    # Define color mapping function for the degree
    col_fun <- colorRamp2(
      range(degree_vals),
      c("#ffe577", "#ff5a00")
    )

    col_fun2 <- colorRamp2(
      range(closeness_vals),
      c("#005f73", "#94d2bd")
    )

    col_fun3 <- colorRamp2(
      range(eigenvector_vals),
      c("#d0ada7", "#5d2e46")
    )

    # Add as row annotation
    row_anno <- rowAnnotation(
      degree = degree_vals,
      closeness = closeness_vals,
      eigenvector = eigenvector_vals,
      col = list(degree = col_fun,
                 closeness = col_fun2,
                 eigenvector = col_fun3)
    )

    # Generate heatmap with cluster annotations
    complex_heatmap = Heatmap(M2,
                              cluster_rows = row_hclust,
                              cluster_columns = col_hclust,
                              show_row_names = TRUE,
                              show_column_names = TRUE,
                              top_annotation = col_anno,
                              left_annotation = row_anno,
                              right_annotation = right_anno,
                              heatmap_legend_param = list(title = "Matrix Values"),
                              row_title = "",
                              column_title = "",

                              cell_fun = function(j, i, x, y, width, height, fill) {
                                if(annotate_matrix[rownames(M2)[i],colnames(M2)[j]]){
                                  grid.points(x, y, pch = 16, size = unit(3, "mm"))
                                }
                              }

    )


  }else{
    # Generate heatmap with cluster annotations
    complex_heatmap = Heatmap(M2,
                              cluster_rows = row_hclust,
                              cluster_columns = col_hclust,
                              show_row_names = TRUE,
                              show_column_names = TRUE,
                              # right_annotation = row_anno,
                              top_annotation = col_anno,
                              heatmap_legend_param = list(title = "Matrix Values"),
                              row_title = "Rows (Samples)",
                              column_title = "Columns (Features)"
                              # col = colorRamp2(c(min(M2), 0, max(M2)), c("blue", "yellow", "red"))
    )
  }

  return(list(#"heatmap_analysis" = heatmap_analysis,
    "clusters" = clusters,
    "n_genes_par_cluster" = n_genes_par_cluster,
    "average_correlation_in_clusters"=average_correlation_in_clusters,
    "list_correlation_in_clusters"=list_correlation_in_clusters,
    "complex_heatmap"=complex_heatmap))
}

#' Create a D3 Network Visualization
#'
#' This function generates a D3 network visualization based on statistical data and a graph structure.
#'
#' @param statistics A data frame containing the statistics on the gene-gene correlation data.
#' @param g An igraph graph object representing gene-gene interactions.
#' @param th A threshold for filtering correlations.
#' @param positive Logical, whether to consider positive correlations.
#'
#' @return A list containing the D3 network visualization, a data frame of node statistics, and an igraph graph object.
#'
#' @export
make_d3_network = function(statistics,g, th = 0, positive = TRUE){
  # correlation = statistics[,"CorGenePatteerns"]
  g_edges = E(g)
  g_edges_names = attr(g_edges,"vnames")
  g_edges_type =E(g)$type

  mat = statistics[,c("Feature 1","Feature 2","CorGenePatteerns")]
  to_rem = which(abs(mat[,3]) < th)
  if (length(to_rem) > 0) mat = mat[-to_rem,]

  ug = unique(as.vector(as.matrix(mat[,1:2])))
  ug = ug[ug %in% V(g)$name]
  if (length(ug) == 0) return(NULL)

  pb = txtProgressBar(min = 0, max = nrow(mat), style = 3)
  Adj_neg = c()
  for (i in 1:nrow(mat)) {

    if (sum(c(mat[i,1], mat[i,2]) %in% V(g)$name) < 2){
      ei = 0
    }else{
      ei <- igraph::get.edge.ids(g, vp = c( mat[i,1], mat[i,2]))
    }
    if (ei != 0) {
      res = strsplit(x = g_edges_names[ei],split = "\\|")
      x = do.call(rbind,res)
      x = cbind(x,g_edges_type[ei])
      x = cbind(x, mat[i,3])
      Adj_neg = rbind(Adj_neg,x)
    }

    setTxtProgressBar(pb,i)
  }
  close(pb)

  if (is.null(Adj_neg)) {
    print("There is no edges in between the pair of genes on the PPI!
          Edges will be drawn based on the correlation patterns")

    mat2 = mat
    mat2 = mat2[mat2[,3]>0,]
    el = as.matrix(mat2[,1:2])
    rownames(el) = NULL
    g_mat = igraph::graph_from_edgelist(el = el)
    E(g_mat)$weight = as.numeric(abs(mat2[,3]))
    E(g_mat)$type = "correlation"

    node_degree = igraph::degree(g_mat)
    node_betwenness = igraph::betweenness(g_mat)
    df = data.frame(degree=node_degree, betwenness = node_betwenness)

    g_mat_d3 = networkD3::igraph_to_networkD3(g_mat)

  }else{
    colnames(Adj_neg) = c("from","to","type","correlation")

    g_mat_pos = igraph::graph_from_edgelist(el = Adj_neg[,1:2])
    E(g_mat_pos)$weight = abs(as.numeric(Adj_neg[,"correlation"]))
    E(g_mat_pos)$type = Adj_neg[,"type"]
    g_mat = igraph::graph_from_edgelist(el = Adj_neg[,1:2])

    E(g_mat)$weight = as.numeric(Adj_neg[,"correlation"])
    E(g_mat)$type = Adj_neg[,"type"]

    node_degree = igraph::degree(g_mat)
    node_betwenness = igraph::betweenness(g_mat_pos)
    df = data.frame(degree=node_degree, betwenness = node_betwenness)

    g_mat_d3 = networkD3::igraph_to_networkD3(g_mat)
  }


  Group = rep("",nrow(g_mat_d3$nodes))
  g_mat_d3$nodes = cbind(g_mat_d3$nodes,Group)
  value = E(g_mat)$weight
  g_mat_d3$links = cbind(g_mat_d3$links,value)

  mycolors = c("green","red","blue","orange")
  names(mycolors) = c("ppi","tfs","regulation","correlation")

  if(nrow(g_mat_d3$links)>0){

    p = networkD3::forceNetwork(Links = g_mat_d3$links,Nodes = g_mat_d3$nodes,
                     NodeID = "name",Group="Group",zoom = T, Value = "value",
                     linkColour = mycolors[E(g_mat)$type], legend=F)

  }else{
    p=NULL
  }

  return(list(p=p, df=df,g_mat=g_mat))
}

#' Filter Edges Between Existing Nodes in a Graph
#'
#' This function filters an edge data frame to include only those edges
#' for which both nodes exist in the provided `igraph` object.
#'
#' @param graph An `igraph` object containing the network structure.
#' @param edge_df A data frame where the first two columns represent
#' source and target node names, respectively.
#'
#' @return A filtered data frame containing only the edges for which
#' both nodes exist in the graph.
#' @export
filter_existing_nodes <- function(graph, edge_df,feature1 = "Feature 1", feature2 = "Feature 2") {
  nodes <- V(graph)$name
  edge_df[edge_df[[feature1]] %in% nodes & edge_df[[feature2]] %in% nodes, ]
}

#' Check Presence of Edges in a Graph
#'
#' This function checks whether edges in a data frame are present in
#' the specified `igraph` object. It uses the `filter_existing_nodes()`
#' helper to ensure only valid nodes are checked, and then marks
#' presence via a logical column.
#'
#' @param graph An `igraph` object representing the graph to check against.
#' @param edge_df A data frame where the first two columns are node names
#' for the edges to be checked.
#'
#' @return A data frame with an additional logical column indicating
#' whether each edge is present in the graph.
#' @export
check_edges_in_graph <- function(graph, edge_df, feature1 = "Feature 1", feature2 = "Feature 2") {

  current_edges = paste(edge_df[,feature1],edge_df[,feature2],sep = "|")
  rownames(edge_df) = current_edges

  edge_df_2 = filter_existing_nodes(graph,edge_df,feature1,feature2)
  current_edges = rownames(edge_df_2)

  edge_keys <- as_ids(E(graph))

  presence = current_edges  %in% edge_keys
  names(presence) = current_edges

  print("reverse")
  res  = unlist(lapply(strsplit(x = names(presence)[presence == T],split = "\\|"), function(elem){paste(elem[2], elem[1],sep ="|")}))
  presence[res] = TRUE

  edge_df = cbind(edge_df,presence[rownames(edge_df)])

  return(edge_df)
}





