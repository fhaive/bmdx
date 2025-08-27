# given a model the function estimate the BMD, BMDL, BMDU and AC50 values

#'
#' Given a model, this function estimates the Benchmark Dose (BMD),
#' Benchmark Dose Lower Confidence Limit (BMDL), Benchmark Dose Upper Confidence
#' Limit (BMDU), and AC50 values.
#'
#' @param model The model object representing the dose-response relationship.
#' @param deviation_type Character string specifying the type of deviation from
#'        the fitted model to use for BMD calculation. Default is "standard".
#' @param rl The relative level used to calculate the BMD. Default is 1.349.
#' @param confidence_interval The confidence level for the confidence interval.
#'        Default is 0.95.
#'
#' @return A modified model object with additional attributes for BMDL, BMDU, and AC50.
#'         If an error occurs during the estimation, NA values are assigned to BMDL and BMDU.
#' @export
point_of_departure = function(model,
                              deviation_type = "standard",
                              rl = 1.349,
                              confidence_interval = 0.95){

  # model = suppressWarnings( BMD(model, deviation_type, rl = rl))
  model = BMD(model, deviation_type, rl = rl)

  x_range = seq(min(model$data_x), max(model$data_x), length.out = 1000)
  newdata = data.frame(x_range)
  colnames(newdata) = model$x

  res = tryCatch({
    # conf_interval = suppressWarnings(predict(model, newdata = newdata, interval = "confidence", level = confidence_interval))
    conf_interval = predict(model, newdata = newdata, interval = "confidence", level = confidence_interval)

    if (model$adverse_direction == -1) {
      bmdl <- stats::approx(x = conf_interval[,2], y = newdata[,1], xout =  model$BMR)$y
      bmdu <- stats::approx(x = conf_interval[,3], y = newdata[,1], xout =  model$BMR)$y
    } else {
      bmdl <- stats::approx(x = conf_interval[,3], y = newdata[,1], xout =  model$BMR)$y
      bmdu <- stats::approx(x = conf_interval[,2], y = newdata[,1], xout =  model$BMR)$y
    }

    c("bmdl" = bmdl, "bmdu" = bmdu)

  }, warning = function(w){
    c("bmdl" = NA, "bmdu" = NA)
  },
  error = function(e) {
    c("bmdl" = NA, "bmdu" = NA)
  })

  model[["BMDL"]] = unname(res["bmdl"])
  model[["BMDU"]] = unname(res["bmdu"])
  model = AC50(model)

  return(model)
}


#' Calculate the benchmark response level
#'
#' This function calculates the benchmark response level based on the given risk factor,
#' whether the response should be increased or decreased, and the background level.
#' The calculation is based on Equation 14 from the paper referenced in the code.
#'
#' @param risk_factor The risk factor used to determine the benchmark response level.
#'        Default is 0.1.
#' @param increase Logical value indicating whether the response should be increased.
#'        If TRUE, the response is increased; if FALSE, it is decreased. Default is TRUE.
#' @param background_level The background level of the response. Default is 0.01.
#'
#' @return The benchmark response level calculated based on the input parameters.
#'
#' @export
bmr_factor = function(risk_factor=0.1, increase=TRUE, background_level = 0.01){
  # Equation 14 from doi:10.1111/j.1539-6924.1995.tb00095.x
  response_level = stats::qnorm(background_level/2, lower.tail = FALSE) - stats::qnorm(risk_factor + background_level, lower.tail = FALSE)
  if (!increase) { response_level = -response_level}
  return(response_level)
}


#' Compute deviations for BMD modeling
#'
#' This function computes deviations for BMD modeling based on the specified deviation type.
#' Three different types of deviation are available: standard deviation, relative deviation,
#' and absolute deviation. Refer to the EPA documentation for detailed information on the
#' deviation types (slide 7) from the provided link.
#'
#' @param deviation_type Character string specifying the type of deviation to compute.
#'        Possible values are "standard", "relative", and "absolute".
#' @param model The model object used for computation.
#' @param rl The relative level used in the deviation calculation. Default is 1.349.
#' @return The computed deviation value based on the specified deviation type.
#'
#' @export


compute_deviations = function(deviation_type = "standard", model, rl = 1.349){#, constant_variance = TRUE,model_variance=FALSE,significance_level=0.05){
  # three different types of deviation are available. Namely the Standard deviation, relative deviation and absolute deviation.
  # Refer to EPA documentation for details: (slide 7) of https://clu-in.org/conf/tio/bmds/slides/BMDS_Continuous_Models.pdf

  starting_point = predict(model, data.frame(dose = 0))

  if (deviation_type == "standard") {
    respLev = rl * sqrt(model$control_variance) #sqrt(control_variance(model, constant_variance = constant_variance, model_variance = model_variance,significance_level = significance_level))
  }
  if (deviation_type == "relative") {
    respLev = rl * starting_point
  }
  if (deviation_type == "absolute") {
    respLev = rl
  }

  diff_resp = starting_point + (respLev * model$adverse_direction)
  return(diff_resp)
}
