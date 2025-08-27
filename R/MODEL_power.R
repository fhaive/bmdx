predict.bmdx.power_model = function(x, ...){ as.matrix(investr::predFit(x$fitted, ...)) }
residuals.bmdx.power_model = function(x){ stats::residuals(x$fitted) }

power_model = function(x,y, g=0, b=1, d=1, max_iter = 250){

  fo = paste(y, "~","g+b*", x,"^d",sep = "")

  param_list = list("x" = x,"y" = y, "params" = c("g" = g,"b" = b,"d" = d),
                    "formula" = fo,
                    mod_name = "power",
                    max_iter = max_iter)

  class(param_list) = c("bmdx.power_model", "bmdx")
  return(param_list)
}


print.bmdx.power_model = function(model){
  utils::str(model)
  invisible(model)
}


fit.bmdx.power_model = function(model, data, weights){

  lower_mod = c(g = -Inf, b = -Inf , d = 0)
  upper_mod = c(g = Inf, b = Inf , d = 18)

  args <- list(
    formula = model$formula,
    data = data,
    start = model$params,
    algorithm = "LM",
    lower = lower_mod,
    upper = upper_mod,
    control = minpack.lm::nls.lm.control(maxiter = model$max_iter),
    model = TRUE
  )
  if (!missing(weights)) { args$weights <- weights }
  m <- do.call(minpack.lm::nlsLM, args)

  model[["log_likelihood"]] = unname(stats::logLik(m)[[1]])
  model[["aic"]] = unname(stats::AIC(m))
  model[["residual_variance"]] = unname((summary(m)$sigma)^2)
  model[["R2"]] = modelr::rsquare(m, data)

  model[["fitted"]] = m

  #model[["lack_of_fit"]] = NA #1 #TODO: decide what to do with lack of fit for non linear models
  coef_mono = stats::coef(m)["b"]
  model[["adverse_direction"]] =  ifelse(coef_mono < 0 , -1, 1)

  model[["data_x"]] = data[, model$x]
  model[["data_y"]] = data[, model$y]

  invisible(model)
}


BMD.bmdx.power_model = function(model,deviation_type = "standard", rl = 1.349){#, constant_variance = TRUE,model_variance=FALSE,significance_level=0.05){
  diff_resp = compute_deviations(deviation_type, model, rl)#, constant_variance,model_variance,significance_level)

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

AC50.bmdx.power_model = function(model){

  ac50 = tryCatch({
    doseRange = seq(from = min(model$data_x),to = max(model$data_x), length.out = 1000)
    new_df = data.frame(doseRange)
    colnames(new_df) = model$x
    fitVal = predict(model, new_df)

    half_resp = max(fitVal) - (abs(max(fitVal) - min(fitVal))/2)
    ac50 <-  stats::approx(x = fitVal, y = doseRange, xout = half_resp)$y

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
