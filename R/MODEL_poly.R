
poly_model = function(x,y,degree=2){

  PolyFitList <- list("Linear" = NULL,
                      "Quadratic" = paste(y, "~", x," + I(",x, "*", x,")"),
                      "Cubic" = paste(y, "~", x, "+ I(",x,"*", x,") + I(",x,"*", x,"*", x,")"),
                      "4thgrade" = paste(y, "~", x, "+ I(",x,"*", x,") + I(",x,"*", x,"*", x,") + I(", x,"*", x,"*", x,"*",x,")"),
                      "5thgrade" = paste(y, "~", x, "+ I(",x,"*", x,") + I(",x,"*", x,"*", x,") + I(", x,"*", x,"*", x,"*",x,") + I(",x,"*", x,"*", x,") + I(", x,"*", x,"*", x,"*",x,"*",x,")"))

  param_list = list("x" = x,"y" = y,"degree" = degree,"formula" = PolyFitList[[degree]], mod_name = paste("poly",degree, sep = ""))

  class(param_list) = c("bmdx.poly_model", "bmdx")

  return(param_list)
}


print.bmdx.poly_model = function(model){
  utils::str(model)
  invisible(model)
}

fit.bmdx.poly_model = function(model, data, weights){

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
  coef_mono = as.numeric(stats::coef(stats::lm(data[, model$y]~data[, model$x]))[2])
  model[["adverse_direction"]] =  ifelse(coef_mono < 0 , -1, 1)

  model[["data_x"]] = data[, model$x]
  model[["data_y"]] = data[, model$y]

  invisible(model)
}


BMD.bmdx.poly_model = function(model,deviation_type = "standard", rl = 1.349){
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

predict.bmdx.poly_model = function(x, ...){ as.matrix(investr::predFit(x$fitted, ...)) }
residuals.bmdx.poly_model = function(x){ stats::residuals(x$fitted) }


AC50.bmdx.poly_model = function(model){

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






