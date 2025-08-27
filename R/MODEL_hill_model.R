
hill_model = function(x,y, g=0, v=1, n=1, k=10, max_iter = 500){

  fo = paste(y, "~","g+((v*", x, "^n)/(k^n+", x,"^n))",sep = "")

  param_list = list("x" = x,"y" = y,
                    "params" = c("g" = g,"v" = v,"n" = n, "k" = k),
                    "mod_name" = "hill",
                    "max_iter" = max_iter,
                    "formula" = fo)

  class(param_list) = c("bmdx.hill_model", "bmdx")
  return(param_list)
}


fit.bmdx.hill_model = function(model, data, weights){

  model_lower = c(g = -Inf, v = -Inf , n = 1, k = 1e-3)
  model_upper = c(g = Inf, v = Inf , n = 18, k = Inf)

  args <- list(
    formula = model$formula,
    data = data,
    start = model$params,
    algorithm = "LM",
    lower = model_lower,
    upper = model_upper,
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
  coef_mono = stats::coef(m)["v"]
  model[["adverse_direction"]] =  ifelse(coef_mono < 0 , -1, 1)

  model[["data_x"]] = data[, model$x]
  model[["data_y"]] = data[, model$y]

  invisible(model)
}

predict.bmdx.hill_model = function(x, ...){ as.matrix(investr::predFit(x$fitted, ...)) }
residuals.bmdx.hill_model = function(x){ stats::residuals(x$fitted) }

BMD.bmdx.hill_model = function(model, deviation_type = "standard", rl = 1.349){#, constant_variance = TRUE,model_variance=FALSE,significance_level=0.05){
  # starting_point = predict(model, data.frame(dose = 0))
  #
  # respLev = rl * sqrt(control_variance(model, constant_variance = constant_variance, model_variance = model_variance,significance_level = significance_level))
  # diff_resp = starting_point + (respLev * model$adverse_direction)

  diff_resp = compute_deviations(deviation_type, model, rl)#, constant_variance,model_variance,significance_level)


  doseRange = seq(min(model$data_x), max(model$data_x), length.out = 1000)
  # fitVal <- predict(model, newdata=data.frame(dose=doseRange))
  # bmd_val <-  stats::approx(x = fitVal, y = doseRange, xout = diff_resp)$y

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


AC50.bmdx.hill_model = function(model){

  ac50 = tryCatch({
    ac50 = as.numeric(stats::coef(model$fitted)["k"])
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
