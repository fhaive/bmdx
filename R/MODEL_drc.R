# library(drc)


filter_drc_predict_warnings <- function(w) {
  if ( any( grepl( "Recycling array of length 1 in array-vector arithmetic is deprecated.", w) ) ) invokeRestart( "muffleWarning" )
}


drc_model = function(x,y, model, data_type = c("continuous","binomial")){

  data_type <- match.arg(data_type)

  param_list = list(
    "x" = x,"y" = y,
    "mod_name" = paste0("drc.", class(model)),
    "drc_spec" = model,
    "data_type" = data_type
  )

  class(param_list) = c("bmdx.drc_model", "bmdx")
  return(param_list)
}

fit.bmdx.drc_model = function(model, data, weights){

  fo = paste(model$y, "~", model$x ,sep = "")
  args <- list(
    formula = fo,
    data = data,
    fct = model$drc_spec,
    type = model$data_type
  )
  if (!missing(weights)) { args$weights <- weights }
  # m <- do.call(drc::drm, args)

  m = drc::drm(formula = args$formula,data = args$data, fct = args$fct, type = args$type)

  model[["log_likelihood"]] = unname(stats::logLik(m)[[1]])
  model[["aic"]] = unname(stats::AIC(m))
  model[["residual_variance"]] = suppressWarnings(unname(summary(m)$resVar))

  predicted = suppressWarnings(predict(m, data))
  residual = predicted - data[,model$y]
  model[["R2"]] = 1 - stats::var(residual, na.rm = TRUE)/stats::var(data[,model$y], na.rm = TRUE)

  model[["fitted"]] = m

  #model[["lack_of_fit"]] = NA #1 #TODO: decide what to do with lack of fit for non linear models
  # fit a linear model to estimate the slope of the data
  coef_mono = stats::coef(stats::lm(fo, data = data))[model$x]
  model[["adverse_direction"]] =  ifelse(coef_mono < 0 , -1, 1)

  model[["data_x"]] = data[, model$x]
  model[["data_y"]] = data[, model$y]

  invisible(model)
}

predict.bmdx.drc_model = function(x, ...){
  X = withCallingHandlers( as.matrix(predict(x$fitted, ...)), warning = filter_drc_predict_warnings )
  idx = grep(pattern = "Prediction",x = colnames(X))
  if (length(idx) > 0) {colnames(X)[idx] = "fit"}
  idx = grep(pattern = "Lower",x = colnames(X))
  if (length(idx) > 0) {colnames(X)[idx] = "lwr"}
  idx = grep(pattern = "Upper",x = colnames(X))
  if (length(idx) > 0) {colnames(X)[idx] = "upr"}
  return(X)
}

residuals.bmdx.drc_model = function(x){
  suppressWarnings( stats::residuals(x$fitted))
}


BMD.bmdx.drc_model = function(model,deviation_type = "standard", rl = 1.349){
  diff_resp = compute_deviations(deviation_type, model, rl)

  doseRange = seq(min(model$data_x), max(model$data_x), length.out = 1000)

  bmd_value = tryCatch({
    fitVal <- predict(model, newdata = data.frame(dose = doseRange))
    bmd_val <-  stats::approx(x = fitVal, y = doseRange, xout = diff_resp)$y
    bmd_val
  }, warning = function(w){
    return(NA)
  },
  error = function(e) {
    return(NA)
  })

  model[["BMR"]] = unname(diff_resp)
  model[["BMD"]] = unname(bmd_value)

  return(model)
}


AC50.bmdx.drc_model = function(model){
  ac50 = tryCatch({
    ac50 = drc::ED(model$fitted, 50, interval = "none", display = FALSE)
    ac50
  }, warning = function(w){
    return(NA)
  },
  error = function(e) {
    return(NA)
  })

  model[["AC50"]] = unname(ac50[[1]])
  return(model)
}
