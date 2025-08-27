# library(minpack.lm)
fit = function(x, data, control_variance="constant", significance_level=0.05){

  .fit <- function(x, ...) { UseMethod("fit",x) }
  check_variance(control_variance)
  # if ("bmdx.model_average" %in% class(x) && control_variance != "constant") {
  #     stop("The average model can only be trained with `control_variance='constant'`")
  # }

  model =  tryCatch({
    model <- .fit(x, data)
    model
  }, warning = function(w){
    return("Warning with model")
  },
  error = function(e) {
    # print(e)
    return("Error with model")
  })

  if (class(model)[1] == "character") return(model)

  if (control_variance == "infer") {
    control_variance = infer_variance_type(model, significance_level)
  }

  res = residuals(model)
  if (control_variance == "constant") {
    # MSE - use residuals of all samples
    ctrl_variance = mean(res^2)
  }
  if (control_variance == "nonconstant") {
    # RMSE - use residuals of all control samples (dose=0)
    ctrl_variance = mean(res[model$data_x == 0]^2)
  }

  model[["control_variance_type"]] = control_variance
  model[["control_variance"]] = ctrl_variance

  if (!'bmdx.model_average' %in% class(model)) {
    model[["lack_of_fit"]] = tryCatch({
      lof = test_lack_of_fit(m = model)
      lof
    }, warning = function(w){
      return(NA)
    },
    error = function(e) {
      return(NA)
    })

    model$is_monotonic = is_strictly_monotonic(model)


  }
  return(model)
}

infer_variance_type <- function(model, significance_level) {
  # the null hypothesis is that variances across the groups are the same.
  # To reject the null hypothesis we look for a p-value<0.05 meaning that at least one variance is different from the others (non constant)
  # If we cannot reject the null hypothesis we assume constant variance
  # TODO: we can add a correction level
  if (car::leveneTest(y = model$data_y,group = as.factor(model$data_x))$`Pr(>F)`[1] < significance_level) {
      return("nonconstant")
  } else {
      return("constant")
  }
}

check_variance = function(variance_type) {
    if (!variance_type %in% c("constant", "nonconstant", "infer")) {
        stop("Invalid control variance type")
    }
}

monotonicity = function(x,...){
  UseMethod("monotonicity",x)
}

BMD = function(x,...){
  UseMethod("BMD",x)
}

AC50 = function(x,...){
  UseMethod("AC50",x)

}

predict.bmdx = function(x, ...){ UseMethod("predict.bmdx",x) }
residuals.bmdx = function(x){ UseMethod("residuals.bmdx",x) }

test_lack_of_fit = function(m){
  if ("bmdx.model_average" %in% class(m)) {
      stop("Estimation of lack-of-fit for averaged models is not supported.")
  }
  if ("bmdx.drc_model" %in% class(m)) {
    lof_test = drc::modelFit(m$fitted)
    pval = lof_test$`p value`[2]
  } else {
    data_df = data.frame(x = m$data_x, y = m$data_y)
    colnames(data_df) = c(m$x, m$y)
    f = stats::as.formula(paste0(m$y,"~as.factor(",m$x,")",sep = ""))
    saturated_model = stats::lm(f, data = data_df)
    suppressWarnings({lof_test = stats::anova(m$fitted, saturated_model)})
    pval = lof_test$`Pr(>F)`[2]
  }
  return(pval)
}
