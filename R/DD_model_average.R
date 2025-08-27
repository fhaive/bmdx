#' Add average models to a list of fitted models
#'
#' This function adds average models to a list of fitted models based on the model averaging approach described
#' in "A brief guide to model selection, multimodel inference and model averaging in behavioural ecology using Akaike’s
#' information criterion" by Matthew R. E. Symonds and Adnan Moussalli. The average model is computed using the provided
#' models and added to the list under the name "average". The average model is fitted using the data from the first model
#' in the list. The function returns the updated list of fitted models.
#'
#'
#' @param fitted_models The list of fitted models.
#' @param variance_type variance type, possible values are: constant, nonconstant, infer.
#'
#' @return The updated list of fitted models with the average model added.
#' @export
add_average_models = function(fitted_models, variance_type="constant"){
  pb = utils::txtProgressBar(min = 0, max = length(fitted_models),style = 3)
  for (i in 1:length(fitted_models)) {
    models = fitted_models[[i]]
    models2 = models

    x = models[[1]]$x
    y = models[[1]]$y

    if (length(models2) > 1) {
      tryCatch({
        avg_model = model_average(x, y, models)
        df = data.frame(models[[1]]$data_x, models[[1]]$data_y)
        colnames(df) = c(x,y)
        avg_model = fit(avg_model, df,control_variance=variance_type)
        pod <- point_of_departure(avg_model)
        models2[["average"]] <- pod
        fitted_models[[names(fitted_models)[i]]] <- models2
      }, error = function(e) {
        message("Error in calculating point of departure for model ", i, ": ", e$message)
        # Continue to the next iteration
      })

    }
    utils::setTxtProgressBar(pb,i)
  }
  close(pb)
  return(fitted_models)
}
