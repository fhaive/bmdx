#' Build Model Matrix
#'
#' This function constructs the design matrix for a linear model based on the provided
#' phenotype data, variable of interest, and optional covariates.
#'
#' @param pd The phenotype data as a data frame.
#' @param intercept The value to be used for the intercept (default: -1).
#' @param var.int The variable of interest used in the design matrix.
#' @param covariates A character vector specifying the optional covariates to include in the design matrix (default: NULL).
#' @param verbose Logical; if TRUE, print the formula used for model.matrix (default: TRUE).
#'
#' @return The design matrix for the linear model.
#'
build.model.matrix <- function(pd, intercept = -1, var.int, covariates = NULL, verbose = TRUE) {
  des <- NULL
  if (!is.null(covariates)) {
    f <- stats::formula(paste("~", intercept, "+pd$", var.int, "+pd$", paste(covariates, collapse = "+pd$"), sep = ""))
    if (verbose) print(f)
    des <- stats::model.matrix(f)
    colnames(des) <- gsub(paste("pd\\$", var.int, sep = ""), "", colnames(des))
    for (i in 1:length(covariates))
      # colnames(des) <- gsub( paste("pd\\$", covariates[i], sep = ""), "", colnames(des))
      colnames(des) <- gsub("pd\\$", "", colnames(des))

  }
  else {
    f <- stats::formula(paste("~", intercept, "+pd$", var.int, sep = ""))
    if (verbose) print(f)
    des <- stats::model.matrix(f)
    colnames(des) <- levels(as.factor(pd[,var.int]))
  }
  des
}

#' Differential Gene Expression Analysis
#'
#' This function performs differential gene expression analysis using limma's
#' linear model with specified contrasts and adjustment method.
#'
#' @param data The gene expression data as a data frame or matrix.
#' @param des The design matrix representing the experimental design.
#' @param contrasts A character vector specifying the contrasts for analysis.
#' @param adjust.method The method for p-value adjustment (default: "none").
#'   Options are "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none".
#'
#' @return A list containing topTable results for each contrast specified.
#'
diff.gene.expr <- function(data, des, contrasts, adjust.method) {
  colnames(des) <- make.names(colnames(des))
  cont <- limma::makeContrasts(contrasts = contrasts, levels = des)

  fit <- limma::lmFit(data, des)
  fit2 <- limma::eBayes(limma::contrasts.fit(fit, cont))

  list.top.tables <- list()
  for (i in 1:length(contrasts)) {
    list.top.tables[[i]] <- limma::topTable(fit2, number = nrow(data),coef = i, sort.by = "logFC",adjust.method = adjust.method)
  }
  names(list.top.tables) <- contrasts

  list.top.tables
}

#' Perform Differential Expression Analysis
#'
#' This function filter the data based on their differential analysis
#'
#' @param data_dictionary A list of data frames representing the data dictionary.
#' @param experimental_data A list of data frames containing the experimental data.
#' @param phTable A list of data frames representing the phenotype table.
#' @param time_point_variable The name of the variable representing time points.
#' @param dose_variable The name of the variable representing doses.
#' @param samples_variable The name of the variable representing samples.
#' @param fcAdjustment The type of fold change adjustment (default: "Nominal").
#'   Options are "none" or "whatever_other_adjustment".
#' @param fcPval.th The p-value threshold for fold change filtering (default: 0.05).
#' @param fc.th The fold change threshold (default: 1.5).
#' @param nCores The number of cores for parallel processing (default: 1).
#' @param x The name of the x-axis variable for plotting (default: "dose").
#' @param y The name of the y-axis variable for plotting (default: "expr").
#' @param other_variables_id_col A vector of other variable names used for filtering.
#'
#' @return A list containing the fold change results (filtered and unfiltered),
#'   the filtered data dictionary, and the filtered fold change results.
#' @export
perform_differential_expression_analysis_filtering = function(data_dictionary,
                      experimental_data,
                      phTable,
                      time_point_variable,
                      dose_variable,
                      samples_variable,
                      fcAdjustment = "Nominal",
                      fcPval.th = 0.05,
                      fc.th = 1.5,
                      nCores = 1,
                      x="dose",
                      y = "expr",
                      other_variables_id_col = NULL){

  if (fcAdjustment == "Nominal") fcAdjustment = "none"

  FC_list = list()

  for (i in 1:length(experimental_data)) {
    gex = experimental_data[[i]]
    pheno = phTable[[i]]

    rel.var = as.numeric(as.vector(pheno[,dose_variable]))
    idx_ctrl = which(rel.var == 0)
    rel.var[idx_ctrl] = "control"
    rel.var[-idx_ctrl] = "treatment"

    pheno = cbind(pheno, rel.var)
    time = unique(pheno[,time_point_variable])

    for (ti in time) {
      ph = pheno[pheno[,time_point_variable] == ti,]
      gi = gex[,ph[,samples_variable]]
      design = build.model.matrix(pd = ph, intercept = -1,
                                  var.int = "rel.var",
                                  covariates = dose_variable, verbose = TRUE)


      deg <- diff.gene.expr(data = gi, des = design,
                            contrasts = "treatment-control",
                            adjust.method = fcAdjustment)

      FC_list[[paste(names(experimental_data)[i],"_",ti,sep = "")]] = deg[[1]]
    }
  }

  splitted_names = strsplit(x = names(data_dictionary),split = "_")
  nelem = length(splitted_names[[1]])

  pb = utils::txtProgressBar(min = 1, max = length(data_dictionary), style = 3)
  fold_change_dataframe = matrix(NA, nrow = length(data_dictionary), ncol = nelem + ncol(FC_list[[1]]))

  for (i in 1:length(data_dictionary)) {
    riga = unlist(c(splitted_names[[i]], FC_list[[paste(splitted_names[[i]][1],splitted_names[[i]][2],sep = "_")]][splitted_names[[i]][nelem],]))
    fold_change_dataframe[i, ] = riga

    utils::setTxtProgressBar(pb,i)
  }
  close(pb)

  colnames(fold_change_dataframe)  = c("Exp","Time",other_variables_id_col, "Feature", "logFC","AveExpr","t",  "P.Value", "adj.P.Val" ,"B")
  fold_change_dataframe = as.data.frame(fold_change_dataframe)
  fold_change_dataframe = cbind(fold_change_dataframe, usedFilteringPval = fold_change_dataframe$adj.P.Val)

  fold_change_dataframe$logFC = as.numeric(as.vector(fold_change_dataframe$logFC))
  fold_change_dataframe$AveExpr = as.numeric(as.vector(fold_change_dataframe$AveExpr))
  fold_change_dataframe$t = as.numeric(as.vector(fold_change_dataframe$t))
  fold_change_dataframe$P.Value = as.numeric(as.vector(fold_change_dataframe$P.Value))
  fold_change_dataframe$adj.P.Val = as.numeric(as.vector(fold_change_dataframe$adj.P.Val))
  fold_change_dataframe$B = as.numeric(as.vector(fold_change_dataframe$B))
  fold_change_dataframe$usedFilteringPval = as.numeric(as.vector(fold_change_dataframe$usedFilteringPval))

  to_rem = which(colnames(fold_change_dataframe) %in% c("AveExpr","t","B"))
  if (length(to_rem) > 0) {
    fold_change_dataframe = fold_change_dataframe[,-to_rem]
  }
  #preform filtering
  good_index = filter_pval_th(fold_change_dataframe, fcPval.th = fcPval.th, fc.th = log2(fc.th))
  pass_filter = 1:nrow(fold_change_dataframe) %in% good_index
  fold_change_dataframe = cbind(fold_change_dataframe, pass_filter)

  filtered_data_dictionary = data_dictionary[good_index]
  fold_change_data_list_filtered  = fold_change_dataframe[good_index,]

  return(list("fc_res_dataframe" = fold_change_data_list_filtered,
              "filtered_data_dictionary" = filtered_data_dictionary,
              "fc_res_unfiltered_dataframe" = fold_change_dataframe))

}

#' Filter Fold Change DataFrame by p-value and fold change threshold
#'
#' This function filters a fold change data frame based on given p-value and fold change thresholds.
#'
#' @param fold_change_dataframe A data frame containing fold change results.
#' @param fcPval.th The p-value threshold for filtering (default: 0.05).
#' @param fc.th The fold change threshold for filtering (default: 0.58).
#'
#' @return A numeric vector containing the indices of rows that pass the filtering criteria.
#' @export
filter_pval_th = function(fold_change_dataframe, fcPval.th = 0.05, fc.th = 0.58) {

  good_index = which(abs(as.numeric(as.vector(fold_change_dataframe[,"logFC"]))) > fc.th & as.numeric(as.vector(fold_change_dataframe[,"adj.P.Val"])) < fcPval.th)

  return(good_index)

}
