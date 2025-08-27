#' Scale Numbers to the Range 0-1
#' This function scales numeric values to the range 0-1 by dividing each value by the maximum value in the input vector.
#'
#' @param x A numeric vector to be scaled.
#' @return A numeric vector with scaled values in the range 0-1.
#' @export
scale_numbers = function(x){
  x = x/max(x)
  return(x)
}

#' Make ECDF or Histogram Plots for BMD Data
#'
#' This function creates ECDF (Empirical Cumulative Distribution Function) or Histogram plots
#' for the BMD (Benchmark Dose) data in the provided model statistics data frame.
#'
#' @param mod_stats The model statistics data frame.
#' @param rel_variable The column name in \code{mod_stats} representing the experimental condition or treatment (default: "Experiment").
#' @param group_by The column name in \code{mod_stats} used to group the data for faceting (default: NULL).
#' @param is_group_by_numeric Logical; if TRUE, treat the \code{group_by} variable as numeric (default: TRUE).
#' @param other_variables A character vector specifying additional variables for faceting (default: NULL).
#' @param number_of_column The number of columns in the facet grid (default: 2).
#' @param scaling Logical; if TRUE, scale the BMD values (default: TRUE).
#' @param filter_column The column name(s) in \code{mod_stats} to use for filtering (default: "Model").
#' @param filter_by A list of character vectors containing filtering criteria for each \code{filter_column} (default: list(c("linear"))).
#' @param plot_type The type of plot to create; either "ecdf" for ECDF or "histogram" for histogram (default: "ecdf").
#' @param linewidth Integer specifing the width of the lines. Default 1.3
#' @param text_size Integer specifing the font size. Default = 14
#'
#' @return A plotly object representing the ECDF or histogram plot.
#' @export
ecdf_plots = function(mod_stats, rel_variable = "Experiment",
                           group_by = NULL,
                           is_group_by_numeric = TRUE,
                           other_variables = NULL,
                           number_of_column = 2,
                           scaling = TRUE, # to use when compounds with different doses are analysed
                           filter_column = c("Model"),
                           filter_by = list(c("linear")),
                           plot_type = "ecdf",
                           linewidth=1.3,
                           text_size = 14){

  if ( !is.null(filter_column)) {
    mod_stats = filter_df(mod_stats, filter_column = filter_column,filter_by = filter_by)
    if (is.null(mod_stats)) {
      print("Filter too stringent! No data to plot", "", type = "info")
      return(make_empty_plot())
    }
  }

  mod_stats2 = cbind(paste(mod_stats[,rel_variable], mod_stats[,group_by],mod_stats[,other_variables],sep = "_"), mod_stats)
  colnames(mod_stats2)[1] = "combExp"
  experiments = unique(as.character(mod_stats2$combExp))

  scaled_data = c()
  gene_list = list()

  for (ee in experiments) {
    bmd_ee = as.numeric(as.vector(mod_stats2[mod_stats2$combExp == ee & is.na(mod_stats2$BMD) == FALSE,"BMD"]))
    genes_ee = mod_stats2[mod_stats2$combExp == ee & is.na(mod_stats2$BMD) == FALSE,"Feature"]

    if (length(bmd_ee) == 0) {
      next
    }

    if (scaling) {
      bmd_ee = scale_numbers(bmd_ee)
    }

    x = unlist(strsplit(x = ee,split = "_"))
    scaled_data = rbind(scaled_data, data.frame(bmd = bmd_ee, gene = genes_ee,experiment = ee, rel_variable = x[1], group_by = x[2], other_variables = x[3]))
  }

  if (is_group_by_numeric) {
    scaled_data$group_by = factor(scaled_data$group_by, levels = sort(unique(as.numeric(scaled_data$group_by )),decreasing = F) )
    if (all(is.na(scaled_data$group_by))) {
      return("Is the group variable numeric?")
    }
  }

  scaled_data <- dplyr::arrange(scaled_data, bmd)
  text = paste0(scaled_data$gene)

  if (plot_type == "ecdf") {
    res2 = ggplot2::ggplot(scaled_data, ggplot2::aes(x = bmd, colour = rel_variable)) +
      ggplot2::stat_ecdf(geom = "step", linewidth = linewidth) +  # Increase ECDF line width
      ggplot2::labs(x = "BMD", y = "ECDF") +
      ggplot2::scale_color_discrete(name = "") +
      ggplot2::theme(text = ggplot2::element_text(size =text_size))  # Increase font size
  } else {
    res2 = ggplot2::ggplot(scaled_data, ggplot2::aes(x = bmd, colour = rel_variable)) +
      ggplot2::geom_histogram(alpha = 0.6, position = 'identity', size = linewidth) +  # Adjust histogram line width
      ggplot2::theme(text = ggplot2::element_text(size = text_size))  # Increase font size
  }

  if (is.null(other_variables)) {
    res2 = res2 + ggplot2::facet_wrap(~group_by, ncol = number_of_column)
  } else {
    res2 = res2 + ggplot2::facet_wrap(~group_by + other_variables, ncol = number_of_column)
  }

  return(res2)

}

#' Create an UpSet Plot
#'
#' This function generates an UpSet plot for the given data, which displays the intersections of sets of genes across different experiments or conditions.
#'
#' @param mod_stats A data frame containing the model statistics and gene information.
#' @param rel_variable The name of the variable representing the experiments or conditions.
#' @param group_by The name of the variable to group the data for generating UpSet plots.
#' @param other_variables Additional variables used for plotting, if any.
#' @param filter_column The name of the column to filter the data.
#' @param filter_by A list of filter values to apply on the specified filter_column.
#' @param nintersects The number of intersections to show in the UpSet plot.
#' @param group.by The variable used to group the intersections (e.g., "degree" or "freq").
#' @param order.by The variable used to order the intersections (e.g., "frequency", or "degree").
#'
#' @return An UpSet plot displaying the intersections of sets of genes across different experiments or conditions.
#' @export
#'
upset_plot = function(
    mod_stats, rel_variable = "Experiment",
    group_by = NULL,
    other_variables = NULL,
    filter_column = c("Model"),
    filter_by = list(c("linear")),
    nintersects = 3,
    group.by = "degree",
    order.by ="degree",
    text.scale = c(1.5,1.5,1.2,1.2,1.5,1.5)){

  if (!is.null(filter_column)) {
    mod_stats = filter_df(mod_stats, filter_column = filter_column,filter_by = filter_by)
    # if(is.null(mod_stats)) return(list("heatmap" = make_empty_plot(),"accumulation_plot"=make_empty_plot(),"upset"=make_empty_plot(), "gVars"=gVars, "logVars"=logVars))

    if (is.null(mod_stats)) {
      return("Filter too stringent! No data to plot")
      return(make_empty_plot())
    }
  }

  # mod_stats2 = cbind(paste(mod_stats[,rel_variable], mod_stats[,group_by],mod_stats[,other_variables],sep = "_"), mod_stats)

  mod_stats2 <- cbind(
    paste0(
      mod_stats[, rel_variable],
      ifelse(mod_stats[, group_by] != "", paste0("_", mod_stats[, group_by]), ""),
      ifelse(mod_stats[, other_variables] != "", paste0("_", mod_stats[, other_variables]), "")
    ),
    mod_stats
  )

  colnames(mod_stats2)[1] = "combExp"
  experiments = unique(as.character(mod_stats2$combExp))

  gene_list = list()

  for (ee in experiments) {
    gene_list[[ee]] = mod_stats2[mod_stats2$combExp == ee & is.na(mod_stats2$BMD) == FALSE,"Feature"]
  }

  if (length(gene_list) == 1) {
    res3 = make_empty_plot()
    return("Not enought data for the upsetplot!")

  }else{
    res3 = UpSetR::upset(UpSetR::fromList(gene_list),
                 nsets = length(gene_list),
                 nintersects = nintersects,
                 group.by = group.by,
                 order.by = order.by,text.scale = text.scale)
  }


  return(res3)
}

#' Compute Gene Frequencies and Create Lollipop Plots
#'
#' This function computes gene frequencies based on the provided model statistics and generates lollipop plots for the top genes based on their frequencies.
#'
#' @param mod_stats A data frame containing the model statistics and gene information.
#' @param th Threshold value for gene percentage. Only genes with a percentage above this threshold will be plotted.
#' @param top_genes Maximum number of genes to plot. Default = 100
#' @param rel_variable The name of the variable representing the experiments or conditions.
#' @param group_by The name of the variable to group the data for generating lollipop plots.
#' @param split_by The name of the variable to split the data and generate separate lollipop plots.
#'
#' @return A list containing gene lists, lollipop plots, and matrices for each split-by value or for all data.
#'
#' @export
compute_gene_frequencies = function(mod_stats,
                                    th = 0.7,
                                    top_genes = 100,
                                    rel_variable = "Experiment",
                                    group_by = "None",
                                    split_by = "None"){
  if (group_by == "None") group_by = NULL
  if (split_by == "None") split_by = NULL

  if (!is.null(split_by)) {
    possible_values = as.character(unique(mod_stats[,split_by]))

    lollipol_plot_list = list()
    line_plot_list = list()

    gene_list  = list()
    m_list = list()

    for (vl in possible_values) {
      if (!is.null(group_by)) {
        mod_stats2 = cbind(paste(mod_stats[mod_stats[,split_by] == vl,rel_variable], mod_stats[mod_stats[,split_by] == vl,group_by],sep = "_"), mod_stats[mod_stats[,split_by] == vl,])
        colnames(mod_stats2)[1] = "combExp"
      }else{
        mod_stats2 = cbind(mod_stats[mod_stats[,split_by] == vl,rel_variable], mod_stats[mod_stats[,split_by] == vl,])
        colnames(mod_stats2)[1] = "combExp"
      }

      experiments = unique(as.character(mod_stats2$combExp))

      if(length(experiments)==1){
        print("no gene frequency can be performed on one experiment only")
        return(NULL)
      }

      M1 = table(mod_stats2[,"combExp"], mod_stats2[,"Feature"])
      variables_of_interests = colnames(M1)
      M = as.matrix(M1)

      M[M > 0] = 1
      index = sort(colSums(M), decreasing = T)
      M = M[,names(index)]
      M = t(M)

      df = data.frame(gene = factor(variables_of_interests, levels = rev(variables_of_interests)),
                      percentage = (rowSums(M)/ncol(M) * 100))

      df = df[df$percentage > th,]

      if (nrow(df) == 0) {
        lollipop_plot = make_empty_plot()
        line_plot = make_empty_plot()
      }else{
        lollipop_plot = ggplot2::ggplot(df, ggplot2::aes(x = gene, y = percentage)) +
          ggplot2::geom_segment(ggplot2::aes(x = gene, xend = gene, y = 0, yend = percentage),color = "gray", lwd = 1) +
          ggplot2::geom_point(size = 4, pch = 21, bg = 4, col = 1) +
          ggplot2::scale_x_discrete(labels = df$gene) +
          ggplot2::coord_flip() +
          ggplot2::theme_minimal()

        line_plot = ggplot(df, aes(x = factor(gene), y = percentage, group = 1,
                                        text = paste("Gene:", gene, "<br>Percentage:", percentage, "%"))) +
          geom_line(color = "blue", size = 1) +  # Line plot
          geom_point(size = 3, color = "red") +  # Highlight points
          labs(title = "Gene Frequency", x = NULL, y = "Percentage") +
          theme_minimal() +
          theme(axis.text.x = element_blank(),  # Remove x-axis labels
                axis.ticks.x = element_blank())  # Remove x-axis ticks
      }

      if (nrow(M) == 1) {
        Mnew = matrix(0, nrow = ncol(M), nrow(M))
        rownames(Mnew) = colnames(M1)
        colnames(Mnew) = rownames(M1)
        Mnew[1:nrow(Mnew), 1:ncol(Mnew)] = M[1:nrow(M), 1:ncol(M)]
      }else{
        Mnew = matrix(0, nrow = nrow(M), ncol(M))
        rownames(Mnew) = rownames(M)
        colnames(Mnew) = colnames(M)
        Mnew[1:nrow(Mnew), 1:ncol(Mnew)] = M[1:nrow(M), 1:ncol(M)]
      }

      freq = rowSums(Mnew)/ncol(Mnew)
      Mnew = cbind(Mnew, freq)

      m_list[[vl]] = Mnew
      lollipol_plot_list[[vl]] = lollipop_plot
      line_plot_list[[vl]] = line_plot

      gene_list[[vl]] = df$gene
    }

  }else{

    if (!is.null(group_by)) {
      mod_stats2 = cbind(paste(mod_stats[,rel_variable], mod_stats[,group_by],sep = "_"), mod_stats)
      colnames(mod_stats2)[1] = "combExp"
    }else{
      mod_stats2 = cbind(mod_stats[,rel_variable], mod_stats)
      colnames(mod_stats2)[1] = "combExp"
    }

    experiments = unique(as.character(mod_stats2$combExp))

    if(length(experiments)==1){
      print("no gene frequency can be performed on one experiment only")
      return(NULL)
    }

    M1 = table(mod_stats2[,"combExp"], mod_stats2[,"Feature"])
    # print(M1)

    M = as.matrix(M1)
    # print(M)

    M[M > 0] = 1
    index = sort(colSums(M), decreasing = T)
    M = M[,names(index)]
    # print(M)

    variables_of_interests = colnames(M)

    M = t(M)
    # print(M)


    percentage = (rowSums(M) / ncol(M) * 100)
    df_all = data.frame(gene = factor(names(percentage), levels = names(percentage)),
                        percentage = percentage)


    df = df_all[df_all$percentage > th,]

    if (nrow(df) > top_genes) {
      df = df[1:top_genes,]
      print("Too many genes to plot. Only the top 100 genes will be printed in the lollipop plot. Complete results are in the table below.")
    }

    if (nrow(df) == 0) {
      lollipop_plot = make_empty_plot()
      line_plot = make_empty_plot()
    }else{
      lollipop_plot = ggplot2::ggplot(df, ggplot2::aes(x = gene, y = percentage)) +
        ggplot2::geom_segment(ggplot2::aes(x = gene, xend = gene, y = 0, yend = percentage), color = "gray", lwd = 1) +
        ggplot2::geom_point(size = 4, pch = 21, bg = 4, col = 1) +
        ggplot2::scale_x_discrete(labels = df$gene) +
        ggplot2::coord_flip() +
        ggplot2::theme_minimal()

      line_plot = ggplot(df, aes(x = factor(gene), y = percentage, group = 1,
                                               text = paste("Gene:", gene, "<br>Percentage:", percentage, "%"))) +
        geom_line(color = "blue", size = 1) +  # Line plot
        geom_point(size = 3, color = "red") +  # Highlight points
        labs(title = "Gene Frequency", x = NULL, y = "Percentage") +
        theme_minimal() +
        theme(axis.text.x = element_blank(),  # Remove x-axis labels
              axis.ticks.x = element_blank())  # Remove x-axis ticks
    }

    if (nrow(M) == 1) {
      Mnew = matrix(0, nrow = ncol(M), nrow(M))
      rownames(Mnew) = colnames(M)
      colnames(Mnew) = rownames(M)
      Mnew[1:nrow(Mnew), 1:ncol(Mnew)] = M[1:nrow(M), 1:ncol(M)]
    }else{
      Mnew = matrix(0, nrow = nrow(M), ncol(M))
      rownames(Mnew) = rownames(M)
      colnames(Mnew) = colnames(M)
      Mnew[1:nrow(Mnew), 1:ncol(Mnew)] = M[1:nrow(M), 1:ncol(M)]
    }

    freq = rowSums(Mnew)/ncol(Mnew)
    Mnew = cbind(Mnew, freq)

    m_list = list("All" = Mnew)
    lollipol_plot_list = list("All" = lollipop_plot)
    line_plot_list = list("All" = line_plot)

    gene_list  = list("All" = df$gene)

  }

  result = list("gene_list" = gene_list,
                "lollipol_plot_list" = lollipol_plot_list,
                "line_plot_list" = line_plot_list,
                "m_list"= m_list)
  return(result)
}
