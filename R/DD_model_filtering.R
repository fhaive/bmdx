#' Filter a list of fitted models based on various criteria
#'
#' This function filters a list of fitted models based on specified criteria, such as
#' thresholds, lack of fit, ratios, missing values, R-squared, and monotonicity.
#' The function returns a filtered list of models that pass the specified criteria.
#'
#' @param fitted_models A list of fitted models to be filtered.
#' @param loofth The threshold for lack of fit. Default is 0.1.
#' @param lower_bound_th The lower bound threshold. Default is 0.1% of the lowest dose.
#' @param upper_bound_th The upper bound threshold. Default is 0.1% of the highest dose.
#' @param bmd_bmdl_th The threshold for the ratio of BMD to BMDL. Default is 20.
#' @param bmdu_bmd_th The threshold for the ratio of BMDU to BMD. Default is 20.
#' @param bmdu_bmdl_th The threshold for the ratio of BMDU to BMDL. Default is 40.
#' @param filter_lower_bound Logical value indicating whether to filter models based on
#'        the lower bound threshold. Default is TRUE.
#' @param filter_upper_bound Logical value indicating whether to filter models based on
#'        the upper bound threshold. Default is TRUE.
#' @param filter_by_lack_of_fit Logical value indicating whether to filter models based on
#'        lack of fit. Default is TRUE.
#' @param ratio_filter Logical value indicating whether to filter models based on ratios.
#'        Default is TRUE.
#' @param bmd_na_filter Logical value indicating whether to filter models with missing
#'        BMD values. Default is TRUE.
#' @param bmdl_na_filter Logical value indicating whether to filter models with missing
#'        BMDL values. Default is TRUE.
#' @param bmdu_na_filter Logical value indicating whether to filter models with missing
#'        BMDU values. Default is TRUE.
#' @param ic50_na_filter Logical value indicating whether to filter models with missing
#'        IC50 values. Default is TRUE.
#' @param r2_filter Logical value indicating whether to filter models based on R-squared.
#'        Default is FALSE.
#' @param r2_th The threshold for R-squared. Default is 0.6.
#' @param filter_by_monotonicity Logical value indicating whether to filter models based on
#'        monotonicity. Default is FALSE.
#' @param filter_by_negative_values Logical value indicating wheather to filter models based on negative estimated effective doses
#' @param filter_by_unordered_values Logival value indicating wheather to filter models for which BMDL<BMD<BMDU is not true
#' @return A filtered list of models that pass the specified criteria.
#'
#' @export
model_filtering = function(fitted_models,
                                loofth = 0.1,
                                lower_bound_th = 0.1,
                                upper_bound_th = 0.1,
                                bmd_bmdl_th = 20,
                                bmdu_bmd_th = 20,
                                bmdu_bmdl_th = 40,
                                filter_lower_bound = TRUE,
                                filter_upper_bound = TRUE,
                                filter_by_lack_of_fit = TRUE,
                                ratio_filter = TRUE,
                                bmd_na_filter = TRUE,
                                bmdl_na_filter = TRUE,
                                bmdu_na_filter = TRUE,
                                ic50_na_filter = TRUE,
                                r2_filter = FALSE,
                                r2_th = 0.6,
                                filter_by_monotonicity = FALSE,
                                filter_by_negative_values = TRUE,
                                filter_by_unordered_values = TRUE
){


  filtered_models = list()
  for (j in 1:length(fitted_models)) {

    models = fitted_models[[j]]
    if (length(models) < 1) {
      next
    }
    good_idx = c()
    for (i in 1:length(models)) {

      mod = models[[i]]

      if (is.null(mod$data_x) || is.null(mod$data_y)) {
        pass = FALSE
      }else{
        pass = inner.check.model(mod,
                                 loofth = loofth,
                                 lower_bound_th = lower_bound_th,
                                 upper_bound_th = upper_bound_th,
                                 bmd_bmdl_th = bmd_bmdl_th,
                                 bmdu_bmd_th = bmdu_bmd_th,
                                 bmdu_bmdl_th = bmdu_bmdl_th,
                                 filter_lower_bound = filter_lower_bound,
                                 filter_upper_bound = filter_upper_bound,
                                 filter_by_lack_of_fit = filter_by_lack_of_fit,
                                 ratio_filter = ratio_filter,
                                 bmd_na_filter = bmd_na_filter,
                                 bmdl_na_filter = bmdl_na_filter,
                                 bmdu_na_filter = bmdu_na_filter,
                                 ic50_na_filter = ic50_na_filter,
                                 # na_filters = na_filters,
                                 r2_filter = r2_filter,
                                 r2_th = r2_th,
                                 filter_by_monotonicity = filter_by_monotonicity,
                                 filter_by_negative_values=filter_by_negative_values,
                                 filter_by_unordered_values=filter_by_unordered_values)
      }

      if (pass)  {
        good_idx = c(good_idx, i)
      }
    }

    good_models = models[good_idx]
    if (length(good_models) == 0) {
      next
    }

    filtered_models[[names(fitted_models)[j]]] = good_models
  }


  return(filtered_models)
}

# return TRUE if model pass the filtering
#' Check if a model passes the filtering criteria
#'
#' This function checks if a model passes the specified filtering criteria. The criteria
#' include thresholds, lack of fit, ratios, missing values, R-squared, and monotonicity.
#' The function returns TRUE if the model passes all the criteria, and FALSE otherwise.
#'
#' @param mod The model to be checked.
#' @param loofth The threshold for lack of fit. Default is 0.1.
#' @param lower_bound_th The lower bound threshold. Default is 0.1% of the lowest dose.
#' @param upper_bound_th The upper bound threshold. Default is 0.1% of the highest dose.
#' @param bmd_bmdl_th The threshold for the ratio of BMD to BMDL. Default is 20.
#' @param bmdu_bmd_th The threshold for the ratio of BMDU to BMD. Default is 20.
#' @param bmdu_bmdl_th The threshold for the ratio of BMDU to BMDL. Default is 40.
#' @param filter_lower_bound Logical value indicating whether to filter the model based on
#'        the lower bound threshold. Default is TRUE.
#' @param filter_upper_bound Logical value indicating whether to filter the model based on
#'        the upper bound threshold. Default is TRUE.
#' @param filter_by_lack_of_fit Logical value indicating whether to filter the model based on
#'        lack of fit. Default is TRUE.
#' @param ratio_filter Logical value indicating whether to filter the model based on ratios.
#'        Default is TRUE.
#' @param bmd_na_filter Logical value indicating whether to filter the model with missing
#'        BMD values. Default is TRUE.
#' @param bmdl_na_filter Logical value indicating whether to filter the model with missing
#'        BMDL values. Default is TRUE.
#' @param bmdu_na_filter Logical value indicating whether to filter the model with missing
#'        BMDU values. Default is TRUE.
#' @param ic50_na_filter Logical value indicating whether to filter the model with missing
#'        IC50 values. Default is TRUE.
#' @param r2_filter Logical value indicating whether to filter the model based on R-squared.
#'        Default is FALSE.
#' @param r2_th The threshold for R-squared. Default is 0.6.
#' @param filter_by_monotonicity Logical value indicating whether to filter the model based on
#'        monotonicity. Default is FALSE.
#'
#' @return TRUE if the model passes the filtering criteria, FALSE otherwise.
#' @export
inner.check.model = function(mod,
                             loofth = 0.1,
                             lower_bound_th = 0.1,
                             upper_bound_th = 0.1,
                             bmd_bmdl_th = 20,
                             bmdu_bmd_th = 20,
                             bmdu_bmdl_th = 40,
                             filter_lower_bound = TRUE,
                             filter_upper_bound = TRUE,
                             filter_by_lack_of_fit = TRUE,
                             ratio_filter = TRUE,
                             bmd_na_filter = TRUE,
                             bmdl_na_filter = TRUE,
                             bmdu_na_filter = TRUE,
                             ic50_na_filter = TRUE,
                             r2_filter = FALSE,
                             r2_th = 0.6,
                             filter_by_monotonicity = FALSE,
                             filter_by_negative_values = TRUE,
                             filter_by_unordered_values = TRUE){

  bmd = mod$BMD
  bmdl = mod$BMDL
  bmdu = mod$BMDU
  ic50 = mod$AC50
  r2 = mod$R2

  fitting.pval = mod$lack_of_fit
  if (is.null(fitting.pval)) {
    fitting.pval = NA
  }

  ## THESE LINES ARE USED TO TEST THE FACT THAT ALL THE NR ARE GREATER THAT ZERO
  ## AND THAT BMDL < BMD < BMDU

  min_dose = 0#min(mod$data_x[mod$data_x > 0])
  max_dose = max(mod$data_x[mod$data_x > 0])
  if (is.na(bmdl)) bmdl_2 = min_dose else bmdl_2 = bmdl
  if (is.na(bmdu)){
    if(is.na(bmd)) bmdu_2 = max_dose else bmdu_2 = bmd
  }else {
    bmdu_2 = bmdu
  }
  if (is.na(bmd)) {
    if(is.na(bmdl)) bmd_2 = min_dose else bmd_2 = bmdl
  }else{
    bmd_2 = bmd
  }

  if (is.na(ic50)) ic50_2 = (max_dose - min_dose)/2 else ic50_2 = ic50
  values_are_ordered = bmd_2 >= bmdl_2 & bmdu_2 >= bmd_2
  #values are not negative
  values_are_not_negative = (bmd_2 >= 0) & (bmdl_2 >= 0) & (bmdu_2 >= 0) & (ic50_2 >= 0)

  pass = TRUE

  if(filter_by_negative_values){
    pass = pass & values_are_not_negative
  }
  # if pass is true it will be an accepted model
  if(filter_by_unordered_values){
    pass = pass & values_are_ordered
  }
  #pass = values_are_not_negative & values_are_ordered

  if (bmd_na_filter) {   # none of the predicted values are NA
    pass = pass & (is.na(bmd) == FALSE)
  }

  if (bmdl_na_filter) {   # none of the predicted values are NA
    pass = pass & (is.na(bmdl) == FALSE)
  }

  if (bmdu_na_filter) {   # none of the predicted values are NA
    pass = pass & (is.na(bmdu) == FALSE)
  }

  if (ic50_na_filter) {   # none of the predicted values are NA
    pass = pass & (is.na(ic50) == FALSE)
  }

  if (filter_by_lack_of_fit) {
    if (is.na(fitting.pval)) { #-9999 is a special number used for non linear models where anova cannot be used to compute lack-of-fit pvalues
      pass_loof_test = TRUE
    }else{
      pass_loof_test = fitting.pval >= loofth
    }
    pass = pass  & pass_loof_test
  }

  if (ratio_filter) {
    pass_ratio_filter = (bmd/bmdl <= bmd_bmdl_th) & (bmdu/bmd <= bmdu_bmd_th) & (bmdu/bmdl <= bmdu_bmdl_th)
    pass = pass & pass_ratio_filter
  }

  if (r2_filter) {
    pass_r2_filter = r2 > r2_th
    pass = pass & pass_r2_filter
  }

  if (filter_by_monotonicity) {
    pass_monotonicity = mod$is_monotonic#is_strictly_monotonic(mod)
    pass = pass & pass_monotonicity
  }

  if (filter_lower_bound) {
    # tolerannce on the bmd minimum and maximum value
    # e.g. remove BMDs are extrapolated lower than the lowest dose (for example BMDs lower than 10% of the lowest dose)
    min_tol = min_dose*lower_bound_th
    bmd_pass_lower_bound = (bmd > min_tol) & (bmdl > min_tol) & (bmdu > min_tol) & (ic50 > min_tol)
    pass = pass & bmd_pass_lower_bound
  }

  if (filter_upper_bound) {
    # e.g. remove BMDs that are extrapolated higher than the highest dose (for example BMDs higher than max_dose - 10% of the max dose)
    bmd_pass_upper_bound = (bmd  < (max_dose - (max_dose*upper_bound_th))) &
      (bmdl < (max_dose - (max_dose*upper_bound_th))) &
      (bmdu < (max_dose - (max_dose*upper_bound_th))) &
      (ic50 < (max_dose - (max_dose*upper_bound_th)))
    pass = pass & bmd_pass_upper_bound

  }
  return(pass)
}


#' Check if a fitted model is strictly monotonic
#'
#' This function checks if a fitted model exhibits strict monotonicity over the range of
#' doses. It returns TRUE if the model is strictly monotonic and FALSE otherwise.
#'
#' @param fittedModel The fitted model to be checked.
#'
#' @return TRUE if the model is strictly monotonic, FALSE otherwise.
#' @export
is_strictly_monotonic = function(fittedModel){

  dose = fittedModel$data_x
  range.length = 1000

  step = (max(dose) - min(dose)) / range.length
  x = seq(min(dose),max(dose),length.out = range.length)

  f1  = predict(fittedModel, newdata = data.frame(dose = x[1:(length(x) - 1)]))
  f2  = predict(fittedModel, newdata = data.frame(dose = x[2:length(x)]))

  deriv = (f2 - f1)/step
  verso = sign(deriv)

  if (all(verso == verso[1])) return(TRUE)
  return(FALSE)

}

