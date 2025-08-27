model_average = function(x,y, models){

  param_list = list(
    "x" = x,"y" = y,
    "mod_name" = "model_average",
    "model_list" = models
  )

  class(param_list) = c("bmdx.model_average", "bmdx")
  return(param_list)
}

residuals.bmdx.model_average = function(model){
  return(model$data_y - predict.bmdx.model_average(model))
}

predict.bmdx.model_average = function(model, ...){
  # REF: A brief guide to model selection, multimodel inference
  # and model averaging in behavioural ecology
  # using Akaike’s information criterion
  # Matthew R. E. Symonds &Adnan Moussalli

  # # average with AIC
  # aic = unlist(lapply(model$fitted, "[[", "aic"))
  # akaike_weights = -(aic - min(aic))/2
  # normalized_weights = as.matrix(exp(akaike_weights) / sum(exp(akaike_weights)))

  normalized_weights = model$"normalized_weights"

  # preds = lapply(X = model$fitted, FUN = predict, ...)

  preds = list()
  for(i in 1:length(model$fitted)){
    preds[[i]] = predict(model$fitted[[i]], ...)
  }

  n_samples = nrow(preds[[1]])
  n_cols = ncol(preds[[1]])
  averaged_pred = array(NA, c(n_samples, n_cols))
  for (i in seq(n_cols)) {
    averaged_pred[, i] = sapply(preds, function(x) {x[, i]}) %*% normalized_weights
  }
  if (ncol(averaged_pred) > 1) {
    colnames(averaged_pred) = c("fit","lwr","upr")
  }

  return(averaged_pred)
}

fit.bmdx.model_average = function(model, data){
  # m = lapply(model$model_list, fit, data = data)
  m = model$model_list

  # average with AIC
  aic = unlist(lapply(m, "[[", "aic"))
  akaike_weights = -(aic - min(aic))/2
  #softmax, larger score to the best model
  normalized_weights = as.matrix(exp(akaike_weights) / sum(exp(akaike_weights)))

  model[["normalized_weights"]] = normalized_weights
  model[["fitted"]] = m
  model[["residual_variance"]] = sum(residuals.bmdx.model_average(model)^2)
  model[["R2"]] = NULL

  coef_mono = as.numeric(stats::coef(stats::lm(data[, model$y]~data[, model$x]))[2])
  model[["adverse_direction"]] = ifelse(coef_mono < 0 , -1, 1)
  model[["data_x"]] = data[, model$x]
  model[["data_y"]] = data[, model$y]

  model$is_monotonic = is_strictly_monotonic(model)

  invisible(model)
}

BMD.bmdx.model_average = function(model, deviation_type = "standard", rl = 1.349){
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

AC50.bmdx.model_average = function(model){
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

  model[["AC50"]] = ac50
  return(model)
}
