
linear_model = function(x,y){
  param_list = list("x" = x,"y" = y, "formula" = paste(y,"~",x), mod_name = "linear")

  class(param_list) = c("bmdx.linear_model", "bmdx")

  return(param_list)
}

residuals.bmdx.linear_model = function(x){ stats::residuals(x$fitted) }
predict.bmdx.linear_model = function(x, ...){ as.matrix(investr::predFit(x$fitted, ...)) }

print.bmdx.linear_model  =  function(model){
  utils::str(model)
  invisible(model)
}


fit.bmdx.linear_model = function(model, data, weights){

  args <- list(
    formula = model$formula,
    data = data,
    weights = NULL
  )
  if (!missing(weights)) { args$weights <- weights }
  m <- do.call(stats::lm, args)

  model[["log_likelihood"]] = unname(stats::logLik(m)[[1]])
  model[["aic"]] = unname(stats::AIC(m))
  model[["residual_variance"]] = unname((summary(m)$sigma)^2)
  model[["R2"]] = summary(m)$r.squared

  model[["fitted"]] = m

  #model[["lack_of_fit"]] = test_lack_of_fit(m)
  coef_mono = as.numeric(stats::coef(m)[2])
  model[["adverse_direction"]] =  ifelse(coef_mono < 0 , -1, 1)

  model[["data_x"]] = data[, model$x]
  model[["data_y"]] = data[, model$y]

  invisible(model)
}


BMD.bmdx.linear_model = function(model, deviation_type = "standard", rl = 1.349){#, constant_variance = TRUE, model_variance = FALSE, significance_level = 0.05){
  diff_resp = compute_deviations(deviation_type, model, rl)#, constant_variance,model_variance,significance_level)

  bmd_value = (diff_resp - model$fitted$coefficients[1])/model$fitted$coefficients[2]

  model[["BMR"]] = unname(diff_resp)
  model[["BMD"]] = unname(bmd_value)
  return(model)
}


AC50.bmdx.linear_model = function(model){

  ac50 = tryCatch({
    doses = model$data_x
    ac50 = max(doses) - (abs(max(doses) - min(doses))/2)
    ac50
  }, warning = function(w){
    return(NA)
  },
  error = function(e) {
    return(NA)
  })


  model[["AC50"]] = unname(ac50)

  return(model)
}






