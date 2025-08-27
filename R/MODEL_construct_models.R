#' Build multiple models based on the given model names.
#'
#' This function constructs a list of model objects based on the specified model
#' names. The available model families include "linear", "hill", "power", "poly2",
#' "poly3", "poly4", "poly5", "exp2", "exp3", "exp4", "exp5", "llog2", "llog3",
#' "llog4", "llog5", "mm2", "weibul12", "weibul13", "weibul14", "weibul22",
#' "weibul23", and "weibul24".
#'
#' @param model_names A character vector specifying the model names to build.
#'                    Default is all available model names.
#' @param max_iter Maximum number of iterations for iterative model fitting.
#'                Default is 1024.
#' @param data_type Type of data, either "continuous" or "binomial".
#'                  Default is "continuous".
#' @param x The predictor variable (independent variable).
#' @param y The response variable (dependent variable).
#'
#' @return model_list_result A list containing the specified model objects.
#' @export
build_models = function(model_names = c("linear","hill","power",
                                           "poly2","poly3","poly4","poly5",
                                           "exp2","exp3","exp4","exp5",
                                           "llog2","llog3","llog4","llog5",
                                           "mm2","weibul12","weibul13","weibul14",
                                           "weibul22","weibul23","weibul24"),
                                            max_iter=1024,
                                            data_type=c("continuous","binomial"), x, y){


  if (data_type == "binomial") {
    if (any(c("linear","hill","power",
         "poly2","poly3","poly4","poly5",
         "exp2","exp3","exp4","exp5") %in% model_names)) {

      print("Only the llog, mm, and weibul model families are allowed with binomial data")
      return(NULL)
    }
  }


  model_list = list(
    "linear"   = linear_model(x = x,y = y),
    "hill"     = hill_model(x = x, y = y, max_iter = max_iter),
    "power"    = power_model(x = x, y = y, max_iter = max_iter),
    "poly2"    = poly_model(x = x,y = y, degree = 2),
    "poly3"    = poly_model(x = x,y = y, degree = 3),
    "poly4"    = poly_model(x = x,y = y, degree = 4),
    "poly5"    = poly_model(x = x,y = y, degree = 5),
    "exp2"     = exp_model(x = x,y = y, n_params = 2, max_iter = max_iter),
    "exp3"     = exp_model(x = x,y = y, n_params = 3, max_iter = max_iter),
    "exp4"     = exp_model(x = x,y = y, n_params = 4, max_iter = max_iter),
    "exp5"     = exp_model(x = x,y = y, n_params = 5, max_iter = max_iter),
    "llog2"    = drc_model(x = x,y = y, drc::LL.2(), data_type = data_type),
    "llog3"    = drc_model(x = x,y = y, drc::LL.3(), data_type = data_type),
    "llog4"    = drc_model(x = x,y = y, drc::LL.4(), data_type = data_type),
    "llog5"    = drc_model(x = x,y = y, drc::LL.5(), data_type = data_type),
    "mm2"      = drc_model(x = x,y = y, drc::MM.2(), data_type = data_type),
    "weibul12" = drc_model(x = x,y = y, drc::W1.2(), data_type = data_type),
    "weibul13" = drc_model(x = x,y = y, drc::W1.3(), data_type = data_type),
    "weibul14" = drc_model(x = x,y = y, drc::W1.4(), data_type = data_type),
    "weibul22" = drc_model(x = x,y = y, drc::W2.2(), data_type = data_type),
    "weibul23" = drc_model(x = x,y = y, drc::W2.3(), data_type = data_type),
    "weibul24" = drc_model(x = x,y = y, drc::W2.4(), data_type = data_type)
  )

  model_list_result = model_list[model_names]
  return(model_list_result)
}
