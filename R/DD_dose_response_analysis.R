#' Perform dose-response analysis on a list of models
#'
#' This function fits a list of models to the same end-point and performs dose-response
#' analysis, including estimation of BMD, BMDL, BMDU, and AC50 values.
#'
#' @param data The data used for fitting the models.
#' @param model_list A list of models to fit to the data.
#' @param deviation_type Character string specifying the type of deviation from the fitted
#'        model to use for BMD calculation. Default is "standard". Allowed values are standard and relative
#' @param rl The relative level used to calculate the BMD. Default is 1.349.
#' @param variance_type Character string specifying the type of variance to use in model
#'        fitting. Default is "constant". Other possible values are "non constant", "model" and "inferred"
#' @param confidence_interval The confidence level for the confidence interval. Default is 0.95.
#' @param significance_level The significance level for model fitting. Default is 0.05.
#'
#' @return A list of models with additional attributes for BMDL, BMDU, and AC50 values.
#'         Models that fail to fit or encounter an error during estimation are excluded
#'         from the returned list.
#' @export
dose_response_analysis = function(
    data,
    model_list,
    deviation_type = "standard",
    rl = 1.349,
    variance_type = "constant",
    confidence_interval = 0.95,
    significance_level = 0.05
){

  for (i in 1:length(model_list)) {

    model_list[[i]] =  tryCatch({
      model = fit(x = model_list[[i]], data,variance_type, significance_level)
      if (class(model)[1] == "character") {
        NA
      }else{
        model = point_of_departure(model,deviation_type, rl, confidence_interval)
        model
      }

    },
    error = function(e) {
      # print(e)
      return(NA)
    })
  }

  toRem = which(is.na(model_list))

  if (length(toRem) > 0) {
    model_list = model_list[-toRem]
  }

  return("models_list" = model_list)
}
