#' Perform ANOVA on data dictionary
#'
#' This function performs ANOVA on a data dictionary, computes p-values, performs adjustment if specified,
#' filters the results based on a significance threshold, and returns the ANOVA results and filtered data
#' dictionary.
#'
#' @param data_dictionary The data dictionary containing the data frames to perform ANOVA on.
#' @param anovaAdjustment The adjustment method for p-values. Options include "Nominal" (no adjustment) and
#'   methods supported by the \code{p.adjust} function. Default is "Nominal".
#' @param anovaPval.th The significance threshold for filtering the ANOVA results. Default is 0.05.
#' @param anovaCores The number of cores to use for parallel computation. Default is 1.
#' @param x The column name for the independent variable. Default is "dose".
#' @param y The column name for the dependent variable. Default is "expr".
#' @param other_variables_id_col Additional column names for other variables to include in the ANOVA.
#' @return anova_res_list A list containing the ANOVA results dataframe, filtered data dictionary, unfiltered ANOVA results dataframe, and a plot list.
#' @export
perform_anova = function(data_dictionary,
                         anovaAdjustment = "Nominal",
                         anovaPval.th = 0.05,
                         anovaCores = 1,
                         x="dose",
                         y = "expr",
                         other_variables_id_col = NULL){

  # plot_list = NULL

  pb = utils::txtProgressBar(min = 1, max = length(data_dictionary), style = 3)
  anova_res = list()
  for (i in 1:length(data_dictionary)) {
    df = data_dictionary[[i]]
    df[,x] = as.factor(df[,x])
    pvalue = unlist(summary(stats::aov(formula = stats::as.formula(paste(y,"~",x,sep = "")), data = df)))["Pr(>F)1"]
    anova_sing_res = c(unlist(strsplit(names(data_dictionary)[i],split = "_")), pvalue)
    anova_res[[i]] = anova_sing_res
    utils::setTxtProgressBar(pb,i)

  }
  close(pb)
  anova_res = do.call(rbind, anova_res)

  coln = c("Exp","Time",other_variables_id_col, "Feature", "Pvalue")
  colnames(anova_res) = coln
  anova_res = as.data.frame(anova_res)
  anova_res$Pvalue = as.numeric(as.vector(anova_res$Pvalue))

  #perform adjustment if needed
  anova_res = pval_adjust(anova_res,  anovaAdjustment)

  #preform filtering
  res_filtering = filter_pval(anova_res,  anovaPval.th)
  anova_res_f =  res_filtering[[1]]

  pass_filter = 1:nrow(anova_res) %in% res_filtering[[2]]
  anova_res = cbind(anova_res, pass_filter)

  filtered_data_dictionary = data_dictionary[res_filtering[[2]]]

  anova_res_list = list("anova_res_dataframe" = anova_res_f,
                        "filtered_data_dictionary" = filtered_data_dictionary,
                        "anova_res_unfiltered_dataframe" = anova_res)
  return(anova_res_list)
              # "anova_plot_list" = plot_list))

}

#' Perform the trend test on multiple datasets.
#'
#' This function performs the trend test on a list of datasets stored in the
#' \code{data_dictionary}. The trend test is performed using the \code{mk.test}
#' function from the \code{trend} package. It calculates the trend p-value for
#' each dataset and adjusts the p-values using the specified \code{trendAdjustment}
#' method. It also performs filtering based on the \code{trendPval.th} threshold.
#'
#' @param data_dictionary A list of datasets, each stored as a data frame in the list.
#' @param trendAdjustment The method for adjusting p-values. Default is "Nominal".
#'                        Available options include "Bonferroni", "Holm", "Hochberg",
#'                        "BH", "BY", "fdr", "none".
#' @param trendPval.th The threshold for filtering p-values. Default is 0.05.
#' @param trendCores The number of CPU cores to use for parallel processing. Default is 1.
#' @param x The column name representing the predictor variable (independent variable) in each dataset.
#' @param y The column name representing the response variable (dependent variable) in each dataset.
#' @param other_variables_id_col The column name containing the identifier for other variables.
#'                                Default is NULL.
#'
#' @return A list containing the trend test results and filtered datasets.
#' @export
perform_trend_test = function(data_dictionary,
                              trendAdjustment = "Nominal",
                              trendPval.th = 0.05,
                              trendCores = 1,
                              x="dose",
                              y = "expr",
                              other_variables_id_col = NULL){

  print("test parallel")
  res = parallel::mclapply(X = 1:length(data_dictionary),FUN = function(i){
    df = data_dictionary[[i]]
    df[[x]] = as.factor(df[[x]])
    df2 = stats::aggregate(df[, 2], list(df[[x]]), stats::median)
    colnames(df2) = c(x,y)
    df2[[x]] = as.factor(df2[[x]])
    trend_test = trend::mk.test(df2[[y]], alternative = "two.sided", continuity = FALSE)
    pvalue = trend_test$p.value
    trend_res = c(unlist(strsplit(names(data_dictionary)[i],split = "_")), pvalue)
    return(trend_res)
  },mc.cores = trendCores)

  print("paralel complete")
  trend_res = do.call(rbind, res)

  coln = c("Exp","Time",other_variables_id_col, "Gene", "Pvalue")
  colnames(trend_res) = coln
  trend_res = as.data.frame(trend_res)
  trend_res$Pvalue = as.numeric(as.vector(trend_res$Pvalue))

  #the functions are in the anova file

  #perform adjustment if needed
  trend_res = pval_adjust(trend_res, trendAdjustment)

  #preform filtering
  res_filtering = filter_pval(trend_res, trendPval.th)
  trend_res_f =  res_filtering[[1]]

  pass_filter = 1:nrow(trend_res) %in% res_filtering[[2]]
  trend_res = cbind(trend_res, pass_filter)

  filtered_data_dictionary = data_dictionary[res_filtering[[2]]]

  trend_res_list = list("trend_res_dataframe" = trend_res_f, "filtered_data_dictionary" = filtered_data_dictionary, "trend_res_unfiltered_dataframe" = trend_res)
  return(trend_res_list)

}



#' Adjust p-values in filtering results
#'
#' This function adjusts the p-values in the filtering results using a specified adjustment method.
#'
#' @param filtering_res The filtering results table.
#' @param adjustment_method The adjustment method for p-values. Options include "Nominal" (no adjustment) and
#'   methods supported by the \code{p.adjust} function. Default is "fdr" (false discovery rate).
#'
#' @return The filtering results table with adjusted p-values and the column "usedFilteringPval" indicating the
#'   p-values used for filtering.
#' @export
pval_adjust = function(filtering_res, adjustment_method = "fdr") {
  if (adjustment_method == "Nominal") {
    print("no adjustment performed")
    filtering_res_adj = cbind(filtering_res, usedFilteringPval = filtering_res$Pval)
  } else {
    adjpvals  = stats::p.adjust(filtering_res$Pvalue, method = adjustment_method)
    filtering_res_adj = cbind(filtering_res, usedFilteringPval = adjpvals)
    filtering_res_adj = cbind(filtering_res_adj, adjustedPval = adjpvals)
  }

  filtering_res_adj$usedFilteringPval = as.numeric(as.vector(filtering_res_adj$usedFilteringPval))
  return(filtering_res_adj)

}

#' Filter results based on p-value threshold
#'
#' This function filters the  results based on a specified p-value threshold. It returns the filtered ANOVA results
#' and the indices of the filtered rows.
#'
#' @param filtering_res The filtering results table.
#' @param Pval.th The p-value threshold for filtering the results. Default is 0.05.
#'
#' @return A list containing the filtered results and the indices of the filtered rows.
#' @export
filter_pval = function(filtering_res, Pval.th = 0.05) {
  good_index = which(filtering_res$usedFilteringPval < Pval.th)
  filtering_res_f = filtering_res[good_index,]
  return(list(filtering_res_f,good_index ))

}

#' Plot a pie chart of filtering results
#'
#' This function generates a pie chart to visualize the filtering results. The pie chart shows the distribution of variables
#' and non-variables based on the filtering p-values in the table. The function takes the table and an
#' p-value threshold as input.
#'
#' @param tab_unfiltered The unfiltered table containing the filtering p-values.
#' @param Pval.th The p-value threshold for determining variable significance.
#' @param title Title of the plot. It is a string. default = ""
#'
#' @return A pie chart visualization of the filtering results.
#' @export
plot_filtering_pie_chart = function(tab_unfiltered,Pval.th, title = ""){
  # x <-  c(sum(as.numeric(tab_unfiltered$usedFilteringPval) < as.numeric(Pval.th)),
  #         sum(as.numeric(tab_unfiltered$usedFilteringPval) >= as.numeric(Pval.th)))
  # labels <-  c("Variable","Non variable")

  x <-  c(sum(as.numeric(tab_unfiltered$pass_filter) == TRUE),
          sum(as.numeric(tab_unfiltered$pass_filter) == FALSE))
  labels <-  c("Significan","Non significant")

  piepercent = round(100*x/sum(x), 2)

  data = data.frame(value = piepercent, labels = labels)
  data = as.data.frame(data)
  p = ggplot2::ggplot(data, ggplot2::aes(x = "", y = value, fill = labels)) +
      ggplot2::geom_bar(stat = "identity", width = 1, color = "white") +
      ggplot2::coord_polar("y", start = 0) +
      ggplot2::geom_text(ggplot2::aes(label = paste(value,"%",sep = "")), position = ggplot2::position_stack(vjust = 0.5)) +
      ggplot2::theme_minimal() +
      ggplot2::ggtitle(title)+
      ggplot2::theme(
        axis.title.x = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        panel.grid = ggplot2::element_blank()
      ) +
      ggplot2::guides(fill = ggplot2::guide_legend(title = NULL))

  p
}


