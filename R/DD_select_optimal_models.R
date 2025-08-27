# use method == Model average, only if average models are already being computed
#' Select optimal models from a list of current models
#'
#' This function selects the optimal models from a list of current models based on the specified method. The supported
#' methods are "AIC" (Akaike Information Criterion) and "Model average". For the "AIC" method, the model with the lowest
#' AIC value is selected for each dataset. For the "Model average" method, if average models are already being computed,
#' the average model is selected; otherwise, the first model in the list is selected. The function returns a list
#' containing the optimal models and a table of computed model statistics.
#'
#' @param current_models The list of current models.
#' @param method The method for selecting the optimal models. Supported values are "AIC" and "Model average".
#' @param time_col_id The identifier for the time column in the model statistics table.
#' @param optional_col_ids Optional identifiers for additional columns in the model statistics table.
#' @param nCores The number of CPU cores to use for parallel computation.
#'
#' @return A list containing the optimal models and the computed model statistics table.
#' @export
select_optimal_models = function(current_models, method = "AIC", time_col_id, optional_col_ids = NULL, nCores = 1){

  if (method == "AIC") {
    optimal_models = list()
    for (i in 1:length(current_models)) {
      opt_idx = which.min(unlist(lapply(current_models[[i]], function(elem) elem$aic)))
      optimal_models[[names(current_models)[i]]] = current_models[[i]][opt_idx]
    }
  }

  if (method == "Model average") {
    optimal_models = list()
    for (i in 1:length(current_models)) {
      if (length(current_models[[i]]) > 1) {
        opt_idx = which(names(current_models[[i]]) == "average")
      }else{
        opt_idx = 1
      }
      if(length(opt_idx)>0){
        optimal_models[[names(current_models)[i]]] = current_models[[i]][opt_idx]
      }
    }
  }

  other_variables_id_col = time_col_id
  BMD_tab_optimal = compute_model_statistics(optimal_models,other_variables_id_col = c(other_variables_id_col, optional_col_ids) ,nCores = nCores)

  return(list("optimal_models" = optimal_models, "BMD_tab_optimal" =  BMD_tab_optimal))
}
