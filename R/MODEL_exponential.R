exp_model = function(x,y, n_params = 2, a = 1, b = 1, c = NA, d = 1, max_iter = 500){
  ModelList <- list("one" = NULL,
                      "two" = paste(y, "~a*exp([sign]b*", x,")", sep = ""),
                      "three" = paste(y, "~a*exp(([sign]b*", x,")^d)", sep = ""),
                      "four" = paste(y, "~a*(c-(c-1) * exp(-b*", x,"))", sep = ""),
                      "five" = paste(y, "~a*(c-(c-1) * exp((-b*", x,")^d))", sep = ""))

  params = c("a" = a, "b" = b)
  if (n_params == 3 || n_params == 5) params = c(params, "d" = d)
  if (n_params == 4  || n_params == 5) params = c(params, "c" = c)

  param_list = list("x" = x,"y" = y,
                    "params" = params,
                    "formula" = ModelList[[n_params]],
                    "n_params" = n_params,
                    "mod_name" = paste("exp",n_params, sep = ""),
                    "max_iter" = max_iter)

  class(param_list) = c("bmdx.exp_model", "bmdx")
  return(param_list)
}

print.bmdx.exp_model = function(model){
  utils::str(model)
  invisible(model)
}


fit.bmdx.exp_model = function(model, data, weights){

  coef_mono = as.numeric(stats::coef(stats::lm(data[, model$y]~data[, model$x]))[2])

  if (is.null(model$c)) {
    c_start = ifelse(coef_mono < 0,0.5,2)
  } else{
    c_start = model$c
  }

  start_list = list("one" = NULL,
                    "two" = list("a" = model$params["a"],"b" = model$params["b"]),
                    "three" = list("a" = model$params["a"],"b" = model$params["b"], "d" = model$params["d"]),
                    "four" = list("a" = model$params["a"],"b" = model$params["b"], "c" = c_start),
                    "five" = list("a" = model$params["a"],"b" = model$params["b"], "c" = c_start, "d" = model$params["d"]))


  if (model$n_params %in% c(2,3)) {
    model$formula = gsub(pattern = "\\[sign\\]",
                         replacement = paste(as.character(sign(coef_mono)),"*",sep = ""),
                         x = model$formula)
  }

  lower_list = list("one" = NULL,
                    "two"   = c("a" = 1e-3,"b" = 1e-3),
                    "three" = c("a" = 1e-3,"b" = 1e-3, "d" = 1),
                    "four"  = c("a" = 1e-3,"b" = 1e-3, "c" = ifelse(coef_mono < 0,1e-3,(1 + 1e-3))),
                    "five"  = c("a" = 1e-3,"b" = 1e-3, "c" = ifelse(coef_mono < 0,1e-3,(1 + 1e-3)), "d" = 1))

  upper_list = list("one" = NULL,
                    "two"   = c("a" = Inf,"b" = Inf),
                    "three" = c("a" = Inf,"b" = Inf, "d" = 18),
                    "four"  = c("a" = Inf,"b" = Inf, "c" = ifelse(coef_mono < 0,(1 - 1e-3),Inf)),
                    "five"  = c("a" = Inf,"b" = Inf, "c" = ifelse(coef_mono < 0,(1 - 1e-3),Inf), "d" = 18))


  args <- list(
    formula = model$formula,
    data = data,
    start = model$params,
    algorithm = "LM",
    lower = lower_list[[model$n_params]],
    upper = upper_list[[model$n_params]],
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
  model[["adverse_direction"]] =  ifelse(coef_mono < 0 , -1, 1)


  model[["data_x"]] = data[, model$x]
  model[["data_y"]] = data[, model$y]

  invisible(model)
}


BMD.bmdx.exp_model = function(model, deviation_type = "standard", rl = 1.349){#, constant_variance = TRUE,model_variance=FALSE,significance_level=0.05){
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


AC50.bmdx.exp_model = function(model){
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

predict.bmdx.exp_model = function(x, ...){ as.matrix(investr::predFit(x$fitted, ...)) }
residuals.bmdx.exp_model = function(x){ stats::residuals(x$fitted) }
