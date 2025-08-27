#' Compute model statistics for fitted models
#'
#' This function computes model statistics for a list of fitted models and returns the results as a data frame.
#'
#' @param fitted_models A list of fitted models.
#' @param other_variables_id_col The name of the column representing the other variables in the model statistics data frame.
#' @param is_parallel Whether to compute the statistics in parallel. Default is TRUE.
#' @param nCores The number of cores to use for parallel computation. Default is 2.
#' @return A data frame containing the computed model statistics.
#' @export

compute_model_statistics = function(fitted_models,
                                         other_variables_id_col,
                                         is_parallel = TRUE,
                                         nCores = 2){

  if (length(fitted_models) == 0) {
    return(NULL)
  }

  if (!is_parallel) {

    modStats = matrix(NA, nrow = sum(unlist(lapply(fitted_models, FUN = length))), ncol = length(unlist(strsplit(split = "_",x = names(fitted_models)[1]))) + 15)
    index = 1
    pb = utils::txtProgressBar(min = 0, max = length(fitted_models), style = 3)
    for (j in 1:length(fitted_models)) {


      # print(paste("j = ", j,sep =""))

      models = fitted_models[[j]]
      if (length(models) < 1) next

      for (i in 1:length(models)) {
        stats_i = get_model_stats(model = models[[i]])
        modStats[index,] = c(unlist(strsplit(split = "_",x = names(fitted_models)[j])), stats_i)
        index = index + 1
        # print(paste("i = ", i,sep =""))
        # print(paste("index = ", index,sep =""))


      }
      utils::setTxtProgressBar(pb,j)
    }
    close(pb)

    to_rem = which(apply(X = modStats, MARGIN = 1, FUN = function(riga) sum(is.na(riga)) == length(riga)))
    if (length(to_rem) > 0) {
      modStats = modStats[-to_rem,]
    }

    colnames(modStats) = c("Experiment",other_variables_id_col,"Feature", names(stats_i))


  }else{
    res = parallel::mclapply( 1:length(fitted_models),FUN = function(j){
      models = fitted_models[[j]]

      if (length(models) > 0) {
        modStats = c()
        for (i in 1:length(models)) {
          stats_i = get_model_stats(models[[i]])
          modStats = rbind(modStats, c(unlist(strsplit(split = "_",x = names(fitted_models)[j])), stats_i))
        }
        return(modStats)
      }

    },mc.cores = nCores)

    modStats = do.call(rbind,res)

    ncol = 2 + length(other_variables_id_col)
    colnames(modStats)[1:ncol] = c("Experiment",other_variables_id_col,"Feature")

  }

  modStats = as.data.frame(modStats)
  modStats$Model = as.factor(modStats$Model)
  modStats$BMR = as.numeric(as.vector(modStats$BMR))
  modStats$BMDL = as.numeric(as.vector(modStats$BMDL))
  modStats$BMD = as.numeric(as.vector(modStats$BMD))
  modStats$BMDU = as.numeric(as.vector(modStats$BMDU))
  modStats$AC50 = as.numeric(as.vector(modStats$AC50))
  modStats$Adverse_direction = as.numeric(as.vector(modStats$Adverse_direction))
  modStats$AIC = as.numeric(as.vector(modStats$AIC))
  modStats$Lack_of_fit = as.numeric(as.vector(modStats$Lack_of_fit))
  modStats$Log_likelihood = as.numeric(as.vector(modStats$Log_likelihood))
  modStats$Residual_variance = as.numeric(as.vector(modStats$Residual_variance))
  modStats$R2 = as.numeric(as.vector(modStats$R2))
  modStats$control_variance = as.numeric(as.vector(modStats$control_variance))
  modStats$control_variance_type = as.factor(modStats$control_variance_type)

  return(modStats)
}

#' given a model the function creates a named vector with all the model statistics
#' @param model an object of class bmdx
#' @return model_stats a named vector with statistics
#' @export
get_model_stats = function(model){
  n_digits = 6
  model_stats = c("Model" = model$mod_name,
                  "BMR" = ifelse(is.null(model$BMR), NA, round(model$BMR,n_digits)),
                  "BMDL" = ifelse(is.null(model$BMDL), NA, round(model$BMDL,n_digits)),
                  "BMD" = ifelse(is.null(model$BMD), NA, round(model$BMD,n_digits)),
                  "BMDU" = ifelse(is.null(model$BMDU), NA, round(model$BMDU,n_digits)),
                  "AC50" = ifelse(is.null(model$AC50), NA, round(model$AC50,n_digits)),
                  "Adverse_direction" = ifelse(is.null(model$adverse_direction),NA, model$adverse_direction),
                  "AIC" = ifelse(is.null(model$aic), NA, round(model$aic,n_digits)),
                  "Lack_of_fit" = ifelse(is.null(model$lack_of_fit), NA, round(model$lack_of_fit,n_digits)),
                  "Log_likelihood" = ifelse(is.null(model$log_likelihood), NA, round(model$log_likelihood,n_digits)),
                  "R2" = ifelse(is.null(model$R2), NA, round(model$R2,n_digits)),
                  "Residual_variance" = ifelse(is.null(model$residual_variance), NA, round(model$residual_variance,n_digits)),
                  "control_variance" = round(unname(model$control_variance),n_digits),
                  "control_variance_type" = model$control_variance_type,
                  "is_monotonic" = model$is_monotonic)
                  # "parameters" = extract_model_parameters(model))
  return(model_stats)
}




extract_model_parameters = function(model){

  model_class = class(model)[1]
  model_name = model$mod_name

  if(model_class == "bmdx.drc_model"){
    params = model$fitted$coefficients
  }

  if(model_name %in% c("linear","poly2","poly3", "poly4","poly5")){
    params = model$fitted$coefficients
  }

  if(model_name %in% c("hill","power","exp2","exp3","exp4","exp5")){
    res = summary(model$fitted)
    params =  res[["cov.unscaled"]][,1]
  }

  collapsed_params = paste(names(params),":",params, sep = "", collapse = ";")
  return(collapsed_params)

}
