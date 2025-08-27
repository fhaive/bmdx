#' Fit dose-response models to data dictionary
#'
#' This function fits dose-response models to a data dictionary using the specified model list. It returns
#' a list of fitted models.
#'
#' @param data_dictionary The data dictionary containing the data frames to fit the models to.
#' @param model_list A list of dose-response models to fit.
#' @param deviation_type The type of deviation to use for BMD calculations. Default is "relative".
#' @param rl The constant value for relative deviation. Default is 1.349.
#' @param confidence_interval The confidence level for the BMD calculation. Default is 0.95.
#' @param variance_type The type of variance to assume for the models. Default is "constant".
#' @param significance_level The significance level for the BMD calculation. Default is 0.05.
#' @param is_parallel Whether to perform the fitting in parallel. Default is FALSE.
#' @param nCores The number of cores to use for parallel computation. Default is 2.
#' @return A list containing the fitted models.
#'
#' @export
fitting_list = function(data_dictionary,
                        model_list,
                        deviation_type="relative",
                        rl = 1.349,
                        confidence_interval = 0.95,
                        variance_type = "constant",
                        significance_level = 0.05,
                        is_parallel = FALSE,
                        nCores = 2){

  # browser()

  if (!is_parallel) {
    # all_fitted_models = list()
    all_fitted_models = vector("list", length(data_dictionary)) # Preallocate the list
    names(all_fitted_models) = names(data_dictionary) # Set names in advance
    pb = utils::txtProgressBar(min = 0, max = length(data_dictionary), style = 3)

    for (i in 1:length(data_dictionary)) {

      df = data_dictionary[[i]]
      bmd_res = dose_response_analysis(data = df,
                                       model_list = model_list,
                                       deviation_type = deviation_type,
                                       rl = rl,
                                       variance_type = variance_type,
                                       confidence_interval = confidence_interval,
                                       significance_level = significance_level)

      all_fitted_models[[names(data_dictionary)[i]]] = bmd_res

      utils::setTxtProgressBar(pb, i)
    }

    return(all_fitted_models)
  }else{
    res = parallel::mclapply(X = 1:length(data_dictionary),FUN = function(i){
      df = data_dictionary[[i]]
      bmd_res = dose_response_analysis(data = df,
                                       model_list = model_list,
                                       deviation_type = deviation_type,
                                       rl = rl,
                                       variance_type = variance_type,
                                       confidence_interval = confidence_interval,
                                       significance_level = significance_level)
      return(list(bmd_res = bmd_res, key = names(data_dictionary)[i]))
    },mc.cores = nCores)

    all_fitted_models_parallel = list()

    for (i in 1:length(res)) {
      all_fitted_models_parallel[[res[[i]]$key]] = res[[i]]$bmd_res
    }
    return(all_fitted_models_parallel)

  }

}
