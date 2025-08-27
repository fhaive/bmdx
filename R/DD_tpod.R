#' Calculate the Percentile of a POD Vector
#'
#' This function computes a specified percentile of a given POD (Point of Departure) vector.
#'
#' @param pod_vector A numeric vector containing POD values.
#' @param percentile A numeric value representing the desired percentile (should be within [0,1]).
#'
#' @return A numeric value representing the calculated percentile of the POD vector.
#' @export
#'

tpod_percentile <- function(pod_vector, percentile){

  # Validate pod_vector to ensure it is numeric
  if (!is.numeric(pod_vector)) {
    stop("The pod_vector should be a numeric vector.")
  }

  # Validate percentile to ensure it's a numeric value within the [0,1] range
  if (!is.numeric(percentile) || any(percentile < 0 | percentile > 1)) {
    stop("The percentile should be a numeric value in the range [0,1].")
  }

  # Calculate the requested percentile value, handling NA values
  percentile_value <- quantile(pod_vector, probs = percentile, na.rm = TRUE)

  # Convert result to numeric
  tPOD <- as.numeric(percentile_value)

  # Return the calculated percentile value
  return(tPOD)
}


#' Calculate the Mean of a POD Vector
#'
#' This function computes the mean of a given POD (Point of Departure) vector, excluding NA values.
#'
#' @param pod_vector A numeric vector containing POD values.
#'
#' @return A numeric value representing the mean of the POD vector.
#' @export


tpod_mean <- function(pod_vector){

  # Validate pod_vector to ensure it is numeric
  if (!is.numeric(pod_vector)) {
    stop("The pod_vector should be a numeric vector.")
  }

  # Calculate the mean value, handling NA values
  mean_value <- mean(pod_vector, na.rm = TRUE)

  # Convert result to numeric
  tPOD <- as.numeric(mean_value)

  # Return the calculated mean value
  return(tPOD)
}



#' Calculate the first mode value using kernel density estimation
#'
#' This function computes the mode of a given vector using kernel density estimation (KDE).
#' The mode is the value that has the highest density in the distribution of the data.
#'
#' @param pod_vector A numeric vector for which the mode is to be calculated.
#'
#' @return A numeric value representing the mode of the input vector.
#'
#' @details
#' The function uses the kernel density estimation technique to estimate the probability density function of the given vector.
#' It then identifies the value that corresponds to the highest density, which is returned as the mode.
#' This approach is useful when the distribution is not unimodal or not normally distributed.
#'
#' @export
tpod_first_mode <- function(pod_vector){

  # Define a function to calculate the mode using kernel density estimation
  # This is a helper function that uses density estimation to find the mode
  mode_kde <- function(pod_vector) {
    d <- density(pod_vector)  # Compute kernel density estimation
    mode_index <- which.max(d$y)  # Find index of maximum density
    d$x[mode_index]  # Return the corresponding value for the mode
  }

  # Calculate the first mode value using kernel density estimation
  mode_value <- mode_kde(pod_vector)

  tPOD = as.numeric(mode_value)  # Ensure the result is returned as numeric

  # Return the calculated mode value as numeric
  return(tPOD)
}


#' Calculate the lowest value of a numeric vector, excluding NA values
#'
#' This function computes the lowest value in the given numeric vector, excluding any `NA` values.
#'
#' @param pod_vector A numeric vector for which the lowest value is to be calculated.
#'
#' @return A numeric value representing the lowest value of the input vector, excluding `NA`s.
#'
#' @details
#' The function uses the `min` function with the `na.rm = TRUE` argument to exclude `NA` values
#' when computing the minimum. It returns the lowest value as a numeric result.
#'
#' @export
tpod_lowest <- function(pod_vector, lowest_method, ratio_threshold, rank){

  if (!is.character(lowest_method) || length(lowest_method) != 1 ||
      !lowest_method %in% c("lowest", "LCRD")) {
    stop("Error: lowest_method must be one of 'lowest', 'LCRD'.")
  }

  pod_vector <- sort(pod_vector)

  # Compute tPOD using the selected method
  tPOD <- switch(lowest_method,
                 "lowest"    = as.numeric(pod_vector[rank]),
                 "LCRD"       = calculate_lcrd(pod_vector,ratio_threshold),
                 stop("Invalid tpod_method specified.")) # Failsafe error handling

  # Return the calculated lowest value as numeric
  return(tPOD)
}

#' Calculate tPOD Using an Accumulation-Based Approach
#'
#' This function calculates the tPOD (Point of Departure) by:
#' \enumerate{
#'   \item Identifying the mode and antimode of the input distribution (via \code{get_mode_antimode}).
#'   \item Filtering the data to values below the antimode.
#'   \item Fitting a shape-constrained additive model (via \code{scam}).
#'   \item Applying the Kneedle algorithm to the smoothed cumulative sum to detect the knee point.
#' }
#'
#' @param pod_vector A numeric vector representing the POD (Point of Departure) values.
#' @param plot Logical. If \code{TRUE}, a plot is generated showing the sorted data, the fitted curve,
#'   a vertical line at the antimode, and another vertical line at the estimated tPOD.
#' @param ... Additional arguments passed to \code{get_mode_antimode}.
#'
#' @return A single numeric value representing the tPOD, identified by the Kneedle method on the fitted cumulative distribution.
#'
#' @details
#' \enumerate{
#'   \item \strong{Sorting}: The function sorts \code{pod_vector} in increasing order.
#'   \item \strong{Accumulation}: It computes the cumulative sum of the sorted values (y) and the log-transformed x-axis (log10 of the sorted values).
#'   \item \strong{Filtering}: Using the estimated \code{antimode} from \code{get_mode_antimode}, data points above the antimode are removed.
#'   \item \strong{Fitting}: A shape-constrained additive model (\code{scam}) is then fitted to the filtered data to smooth the cumulative sums.
#'   \item \strong{Knee Detection}: The Kneedle algorithm (\code{kneedle}) is applied to the fitted curve to detect the knee, which is reported as the tPOD.
#'   \item \strong{Optional Plot}: If \code{plot = TRUE}, the sorted data, fitted curve, and vertical lines marking the antimode and tPOD will be displayed.
#' }
#'
#'
#' @export
tpod_accumulation <- function(pod_vector, plot=FALSE, bw_adjust=3){

    # 1. estimate mode and antimode
    pdf_info <- get_mode_antimode(pod_vector, bw_adjust=bw_adjust)

    # 2. compute accumulation values
    pod_vector_sorted <- sort(pod_vector, decreasing = FALSE)
    y <- cumsum(pod_vector_sorted)
    x <- log10(pod_vector_sorted)

    # 3. filter everything after the antimode
    df <- data.frame(x=x, y=y)
    df <- df[pod_vector_sorted < pdf_info$antimode,]

    # 4. smooth the data
    mod <- scam(
        y ~ s(x, bs = "mpi", k = max(5, round(length(df) / 3))),
        data    = df,
        family  = gaussian()
    )

    # 4. Create fitted values, e.g. for plotting
    df$y_pred <- predict(mod, newdata=df, type="response")

    # 5. Estimate the knee point
    k <- kneedle(df$y_pred, df$x)
    tpod <- 10^k$y
    if(plot){
        dev.new()
        plot(x, y, pch = 19, cex = 0.5,
             xlab = "log10(BMD)", ylab = "Accumulation")
        lines(
              df$x,
              df$y_pred,
              col = "red", lwd = 2
        )
        abline(v=log10(pdf_info$antimode), col="red")
        abline(v=log10(tpod), col="green")
    }
    return(tpod)
}

#' Identify the knee point in two vectors using a simple Kneedle-like approach
#'
#' This function normalizes the input vectors \code{x} and \code{y} to the range \[0, 1\]
#' and calculates the difference \code{(y_norm - x_norm)} at each point. The maximum
#' of this difference is used as the index of the knee point, and the corresponding
#' values of \code{x} and \code{y} at this index are returned.
#'
#' @param x A numeric vector representing the x-coordinates.
#' @param y A numeric vector representing the y-coordinates.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{\code{index}}{The index of the knee point in the original vectors.}
#'   \item{\code{x}}{The x-value at the knee point.}
#'   \item{\code{y}}{The y-value at the knee point.}
#' }
#'
#' @details
#' The vectors \code{x} and \code{y} are first normalized to the \[0, 1\] range.
#' The function then computes the difference \code{(y_norm - x_norm)} at each
#' corresponding index and identifies the knee point as the location where this
#' difference is maximized.
#'
#'
#' @export
kneedle <- function(x, y) {
    x_norm <- (x - min(x)) / (max(x) - min(x))
    y_norm <- (y - min(y)) / (max(y) - min(y))
    curve_diff <- y_norm - x_norm
    knee_index <- which.max(curve_diff)
    return(list(index = knee_index, x = x[knee_index], y = y[knee_index]))
}

#' Identify the first two local maxima and the minimum ("valley") between them in a density estimate
#'
#' This function estimates the density of a numeric vector \code{x}, identifies the first two local maxima
#' (peaks), and locates the lowest point (valley) between those peaks. It returns the x-values for the
#' first peak (mode), the valley (antimode), and a function for interpolating the estimated PDF.
#'
#' @param x A numeric vector from which to estimate a density.
#' @param bw_adjust A numeric multiplier for the bandwidth. The bandwidth is chosen via the "SJ" method and
#'   then multiplied by \code{bw_adjust}.
#'
#' @return A list with elements:
#' \describe{
#'   \item{mode}{The x-value of the first local maximum (peak).}
#'   \item{antimode}{The x-value of the local minimum (valley) between the first two peaks.}
#'   \item{pdf}{A function that interpolates (and extrapolates) the estimated density at arbitrary points.}
#' }
#'
#'
#' @export
get_mode_antimode <- function(x, bw_adjust=3){
    # 1. Estimate the density
    d <- density(x, bw="SJ", adjust=bw_adjust)

    # 2. Find local maxima
    #    A point i is a local maximum if y[i] > y[i-1] and y[i] > y[i+1].
    #    We skip the extreme edges for safety (i in [2, length(y)-1]).
    is.peak <- function(y) {
        n <- length(y)
        pks <- rep(FALSE, n)
        for (i in 2:(n-1)) {
            if (y[i] > y[i-1] && y[i] > y[i+1]) {
                pks[i] <- TRUE
            }
        }
        pks
    }

    peaks_idx <- which(is.peak(d$y))

    if(length(peaks_idx) < 2) {
        stop("Less than two peaks found. Check the data and adjust the bandwitdth")
    }
    peak1_idx <- peaks_idx[1]
    peak2_idx <- peaks_idx[2]

    # 3. Find the valley between those two peaks.
    #    One straightforward way is to find the index of the smallest
    #    density between the two peak indices.
    valley_range     <- d$y[peak1_idx:peak2_idx]
    valley_rel_index <- which.min(valley_range)         # Index within that subset
    valley_idx       <- peak1_idx + valley_rel_index - 1

    # Retrieve the x-values for peaks and valley
    peak1   <- d$x[peak1_idx]
    peak2   <- d$x[peak2_idx]
    valley  <- d$x[valley_idx]

    # 4. Create a density function that interpolates the estimated density
    #    at arbitrary points. The 'rule=2' in approx ensures it uses the
    #    boundary values outside the original range (constant extrapolation).
    densfun <- function(x0) {
        approx(d$x, d$y, xout = x0, rule = 2)$y
    }

    # Return them in whatever format you prefer
    return(list(
        mode=peak1,
        antimode=valley,
        pdf=densfun
    ))
}

#' Generate an accumulation plot with multiple tPOD methods
#'
#' This function generates a cumulative sum plot of the given POD vector, and includes
#' vertical lines indicating the tPODs (points of Departure) for multiple methods.
#' The plot can display multiple tPOD values if they all have the same length.
#' The x-axis can optionally be transformed to a logarithmic scale.
#'
#' @param pod_vector A numeric vector representing the POD (Point of Departure) values.
#' @param tpod_method A character vector containing the names of the tPOD methods used
#'   to calculate the points of Departure.
#' @param pod_value A string that represents the value used to calculate the tPODs,
#'   for example, "BMD" or "BMDL" or "BMDU".
#' @param tPOD A numeric vector containing the points of Departure (tPODs)
#'   corresponding to each method. The length of this vector must match the length of
#'   the `tpod_method` vector.
#' @param xlog A logical value indicating whether to apply a logarithmic transformation
#'   to the x-axis (default is `FALSE`).
#'
#' @return A ggplot object representing the cumulative sum plot with vertical lines for tPODs.
#'
#' @details
#' The function generates an accumulation plot of the sorted POD vector and adds vertical
#' lines at the locations of the provided tPOD values. The plot can display multiple tPOD
#' values (with different methods) if the length of `tPOD` and `tpod_method` are equal.
#' The x-axis can optionally be displayed in logarithmic scale.
#'
#' @import ggplot2
#' @import dplyr
#' @export
tpod_plot <- function(pod_vector, tpod_method, pod_value, tPOD, xlog = FALSE, subtitle = "") {

  # Validate that tPOD and tpod_method have the same length
  if (length(tPOD) != length(tpod_method)) {
    stop("The lengths of tPOD and tpod_method must be the same.")
  }

  # Load required libraries
  library(ggplot2)
  library(dplyr)

  # Prepare data for plotting
  sorted_pod_vector <- sort(pod_vector, decreasing = FALSE)
  df <- data.frame(
    x = if (xlog) log10(sorted_pod_vector) else sorted_pod_vector,
    y = cumsum(sorted_pod_vector)
  )

  # Data frame for vertical lines representing tPOD values
  tPOD_df <- data.frame(
    x = if (xlog) log10(tPOD) else tPOD,
    method = tpod_method
  )

  # Generate the plot
  p <- ggplot(df, aes(x = x, y = y)) +
    geom_line(color = "blue", size = 1) +  # Main accumulation curve
    geom_vline(data = tPOD_df, aes(xintercept = x, color = method), linetype = "solid", size = 1)+  # Vertical lines
    geom_text(data = tPOD_df, aes(x = x, y = max(df$y) * 0.9,
                                  label = round(if (xlog) 10^x else x, 2), color = method),
              angle = 90, vjust = -0.5, hjust = 1, size = 5)

  p <- p +
    theme_minimal() +
    labs(
      title = paste0(subtitle),
      # title = "Accumulation Plot",
      # subtitle =paste0(subtitle),
      x = ifelse(xlog, "Doses", "Doses"),
      y = "Cumulative sum",
      color = paste0("tPOD Method (", pod_value, ")")
    ) +
    theme(
      legend.position = "right",
      legend.box = "horizontal",
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 14),
      axis.title.x = element_text(size = 16),
      axis.title.y = element_text(size = 16),
      plot.title = element_text(size = 18, hjust = 0.5)
    )

  p <- p + scale_x_continuous(
    labels = function(x) round(if (xlog) 10^x else x, 2)  #Not log transformed scale on x axis
  )

  return(p)
}

#' Compute tPOD using various methods and plot the result
#'
#' This function computes the tPOD (Point of Departure) using different methods
#' such as percentile, mean, first mode, lowest, or accumulation After calculating the tPOD,
#' it generates a plot displaying the resulting tPOD value on the accumulation POD distribution.
#'
#' @param model_stats A data frame or list containing the model statistics, including the
#'   POD values (e.g., BMD, BMDL, or BMDU) from which tPOD will be computed.
#' @param pod_value A string indicating which POD value to use from the model statistics.
#'   Possible values are "BMD", "BMDL", and "BMDU" (default is "BMD").
#' @param tpod_method A string indicating the method to compute tPOD. Options include:
#'   "percentile", "mean", "first_mode", "lowest", or "accumulation".
#' @param percentile A numeric value between 0 and 1 representing the desired percentile
#'   (used if `tpod_method` is "percentile").
#'
#' @return A numeric value representing the computed tPOD, according to the selected method.
#'
#' @details
#' The function computes the tPOD using one of the following methods:
#' - **percentile**: Computes tPOD at the specified percentile of the POD values.
#' - **mean**: Computes the tPOD as the mean value of the POD values.
#' - **first_mode**: Computes the tPOD using the first mode (most frequent value) of the POD distribution.
#' - **lowest**: Computes the tPOD as the lowest value of the POD distribution.
#' - **accumulation**: Computes the tPOD using the accumulation method.
#'
#' A plot is also generated showing the accumulation plot of the POD distribution with a vertical line indicating
#' the tPOD value.
#'
#' @import ggplot2
#' @import dplyr
#' @export
tpod_computation <- function(model_stats,
                             tpod_method,
                             percentile,
                             pod_value = "BMD",
                             lowest_method = "lowest",
                             bw_adjust=3,
                             ratio_threshold = 1.66,
                             plot =FALSE,
                             robust = F,
                             rank = 1){

  # Validate tpod_method
  if (!is.character(tpod_method) || length(tpod_method) != 1 ||
      !tpod_method %in% c("percentile", "mean", "first_mode", "lowest", "accumulation")) {
    stop("Error: tpod_method must be one of 'percentile', 'mean', 'first_mode', 'lowest', or 'accumulation'.")
  }

  # Ensure required parameters are provided based on tpod_method
  if (tpod_method == "percentile") {
    if (is.null(percentile) || !is.numeric(percentile) || percentile < 0 || percentile > 1) {
      stop("Error: When tpod_method is 'percentile', percentile must be a numeric value between 0 and 1.")
    }
  }

  # Validate pod_value to ensure it's one of the accepted values
  if (!is.character(pod_value) || sum(pod_value %in%  c("BMD","BMDL","BMDU")) != length(pod_value)) {
    warning("The pod value needs to be 'BMD' or 'BMDU' or 'BMDL'")
  }

  pod_vector = model_stats[[pod_value]]

  #Remove small and large numbers to avoid extreme outliers
  if(robust){
    outliers <- boxplot(pod_vector)
    stats <- outliers$stats
    pod_vector <- pod_vector[(pod_vector <= stats[5]) & (pod_vector >= stats[1])]
  }

  # Compute tPOD using the selected method
  tPOD <- switch(tpod_method,
                 "percentile"    = tpod_percentile(pod_vector, percentile),
                 "mean"          = tpod_mean(pod_vector),
                 "first_mode"    = tpod_first_mode(pod_vector),
                 "lowest"        = tpod_lowest(pod_vector, lowest_method, ratio_threshold, rank),
                 "accumulation"  = tpod_accumulation(pod_vector, bw_adjust,plot= plot),
                 stop("Invalid tpod_method specified.")) # Failsafe error handling

  return(tPOD)
}


#' Apply tPOD Methods to Model Statistics
#'
#' This function computes the tPOD (toxicological point of departure) using various methods
#' and returns a data frame with the computed values. Errors during computation are handled gracefully,
#' with NA returned for failed calculations.
#'
#' @param model_stats A data frame containing model statistics used for tPOD computation.
#' @param tpod_methods_list Character vector; methods to apply, options include "percentile", "mean", "first_mode", "lowest", "accumulation".
#' @param pod_value Character; type of point of departure to use: "BMD", "BMDL", or "BMDU". Default is "BMD".
#' @param percentile Numeric; a value between 0 and 1 indicating the target percentile. Default is 0.95.
#' @param lowest_method Character; method for selecting the lowest value: "lowest" or "LCRD". Default is "lowest".
#' @param bw_adjust Numeric; bandwidth adjustment for density estimation methods. Default is 3.
#' @param robust Logical; whether to use robust statistics where applicable. Default is FALSE.
#' @param plot Logical; whether to generate plots during tPOD computation. Default is FALSE.
#' @param ratio_threshold Numeric; threshold for ratio-based decisions. Default is 1.66.
#' @param rank Numeric; rank selection for the tPOD if multiple candidates exist. Default is 1.
#'
#' @return A data frame containing the computed tPOD values with columns: "Method", "Parameter", and "tPOD".
#'
#' @export

apply_tpod_methods = function(model_stats,
                              tpod_methods_list = c("percentile", "mean", "first_mode", "lowest", "accumulation"),
                              pod_value = "BMD",
                              percentile = 0.95,
                              lowest_method = "lowest",
                              bw_adjust = 3,
                              robust = F,
                              plot = F,
                              ratio_threshold = 1.66,
                              rank = 1){


  parameter_list = c(percentile, "NA", "NA", lowest_method, "NA")
  names(parameter_list) = c("percentile", "mean", "first_mode", "lowest", "accumulation")

  parameter_list = parameter_list[tpod_methods_list]

  statistics = c()
  for(i in 1:length(tpod_methods_list)){

    res <- tryCatch(
      {
        tpod_computation(
          model_stats = model_stats,
          pod_value = pod_value,
          tpod_method = tpod_methods_list[i],
          percentile = percentile,
          bw_adjust = bw_adjust,
          lowest_method = lowest_method,
          ratio_threshold = ratio_threshold,
          robust = robust,
          plot = plot,
          rank = rank
        )
      },
      error = function(e) {
        NA
      }
    )

    statistics = rbind(statistics, c(tpod_methods_list[i],parameter_list[i], res))
  }

  statistics = as.data.frame(statistics)
  colnames(statistics) = c("Method","Parameter","tPOD")
  rownames(statistics) = statistics$Method
  return(statistics)
}

#' Filter Data Based on Antimode of BMD Distribution
#'
#' This function filters `optimal_models_stats` based on the **antimode** of the BMD distribution.
#' It splits the data into two groups:
#' - **Before Antimode**: BMD values **less than or equal to** the antimode.
#' - **After Antimode**: BMD values **greater than** the antimode but **below the 95th percentile**.
#'
#' @param optimal_models_stats A dataframe containing the BMD values and corresponding timepoints.
#' @param timepoints_col_variable A character vector specifying the timepoints column variable to analyze.
#'
#' @return A **list** with two elements:
#' - **$before**: Dataframe containing BMD values **before** the antimode.
#' - **$after**: Dataframe containing BMD values **after** the antimode but below the 95th percentile.
#'
#'
#' @export
filter_antimode <- function(optimal_models_stats, timepoints_col_variable) {
  before_antimode <- list()
  after_antimode <- list()


  # Ensure required parameters are a column in the dataframe
  if (!is.character(timepoints_col_variable) || !(timepoints_col_variable %in%  colnames(optimal_models_stats))) {
    stop("Error: The column variable that you select is not a column in the optimal models stats dataframe")
    }

  timepoints <- unique(optimal_models_stats[[timepoints_col_variable]])


  for (time in timepoints) {
    # Filter data for the current timepoint
    opt_time <- optimal_models_stats[optimal_models_stats[[timepoints_col_variable]] == time, ]
    pod_vector <- opt_time$BMD

    # Get antimode information
    pdf_info <- get_mode_antimode(pod_vector)
    valley <- pdf_info$antimode

    # Filter before and after antimode
    before_antimode[[time]] <- opt_time[opt_time$BMD <= valley, ]
    thr95 <- quantile(pod_vector, 0.95, na.rm = TRUE)
    after_antimode[[time]] <- opt_time[opt_time$BMD > valley & opt_time$BMD < thr95, ]
  }

  # Combine lists into dataframes
  before_antimode <- do.call(rbind, before_antimode)
  after_antimode <- do.call(rbind, after_antimode)

  return(list(before = before_antimode, after = after_antimode))
}

#' Filter Data Based on tpod accumulation for Given Timepoints
#'
#' This function filters data based on the BMD values for specified timepoints. It applies the `tpod_accumulation` function
#' to determine the antimode and then divides the data into two groups:
#' - `before`: data where BMD values are less than or equal to the antimode.
#' - `after`: data where BMD values are greater than the antimode but less than the 95th percentile.
#'
#' @param optimal_models_stats A data frame containing the columns `timepoint` and `BMD`.
#'     This data frame holds the information about the optimal models for different timepoints.
#' @param timepoints_col_variable A character vector specifying the timepoints column variable to analyze.
#' @param plot A logical value indicating whether to generate a plot of the accumulation POD
#'   distribution with a vertical line showing the computed tPOD. Default is `FALSE`.
#'
#' @return A list containing two data frames:
#'     - `before`: data frame with entries where BMD values are less than or equal to the antimode.
#'     - `after`: data frame with entries where BMD values are greater than the antimode but less than the 95th percentile.
#'
#' @importFrom stats quantile
#'
#' @export
filter_tpod_acc <- function(optimal_models_stats, timepoints_col_variable) {
  # Ensure required parameters are a column in the dataframe
  if (!is.character(timepoints_col_variable) || !(timepoints_col_variable %in%  colnames(optimal_models_stats))) {
    stop("Error: The column variable that you select is not a column in the optimal models stats dataframe")
  }

  before <- list()
  after <- list()

  timepoints <- unique(optimal_models_stats[[timepoints_col_variable]])

  for (time in timepoints) {
    # Filter data for the current timepoint
    opt_time <- optimal_models_stats[optimal_models_stats[[timepoints_col_variable]] == time, ]
    pod_vector <- opt_time$BMD

    # Get antimode information
    tpod <- tpod_accumulation(pod_vector, bw_adjust = 3)

    # Filter before and after antimode
    before[[time]] <- opt_time[opt_time$BMD <= tpod, ]
    thr95 <- quantile(pod_vector, 0.95, na.rm = TRUE)
    after[[time]] <- opt_time[opt_time$BMD > tpod & opt_time$BMD < thr95, ]
  }

  # Combine lists into dataframes
  before <- do.call(rbind, before)
  after <- do.call(rbind, after)

  return(list(before = before, after = after))
}

#' Calculate the Lowest Consistent Response Dose (LCRD)
#'
#' This function identifies the lowest consistent response dose (LCRD) from a numeric vector of BMC values.
#' The BMCs are ranked in ascending order, and a consistency check is applied such that the ratio
#' between consecutive BMCs (BMC[n+1] / BMC[n]) must be < 1.66. The first BMC for which all subsequent
#' ratios remain below this threshold is declared the LCRD.
#'
#' @param pod_vector A numeric vector of BMC values.
#' @param ratio_threshold The maximum allowed ratio between consecutive BMCs (default: 1.66).
#' @return A list containing the LCRD, the CRGB (consistent response group of BMCs), and the ranked BMCs.

calculate_lcrd <- function(pod_vector, ratio_threshold) {
  if (!is.numeric(pod_vector)) {
    stop("Input must be a numeric vector of BMC values.")
  }

  # 1. Order the BMC values
  bmc_sorted <- sort(pod_vector)

  # 2. Calculate ratios between consecutive BMCs
  ratios <- bmc_sorted[-1] / bmc_sorted[-length(bmc_sorted)]

  # 3. Find the lowest BMC where all subsequent ratios are < threshold
  for (i in seq_along(ratios)) {
    remaining_ratios <- ratios[i:length(ratios)]
    if (all(remaining_ratios < ratio_threshold)) {
      crgb <- bmc_sorted[i:length(bmc_sorted)]
      lcrd <- crgb[1]
      res_LCRD <- list(LCRD = lcrd, CRGB = crgb, ranked_BMCs = bmc_sorted)
      return(res_LCRD$LCRD)
    }
  }

  warning("No consistent response group found with given threshold.")
  return(NULL)
}


#' Plot BMD Density with Optional tPOD Reference Lines
#'
#' This function generates a density plot of Benchmark Dose (BMD) values and overlays
#' vertical dashed lines for any provided (non-NA) tPOD metrics. The supported metrics
#' include \code{lowest}, \code{percentile}, \code{first_mode}, \code{accumulation}, and \code{mean_value}.
#' A customizable x-axis range can also be specified.
#'
#' @param BMD_values A numeric vector of BMD values.
#' @param lowest Numeric value for the "Lowest" tPOD (optional; default is \code{NA}).
#' @param percentile Numeric value for the "Percentile" tPOD (optional; default is \code{NA}).
#' @param mean_value Numeric value for the "Mean" tPOD (optional; default is \code{NA}).
#' @param accumulation Numeric value for the "Accumulation" tPOD (optional; default is \code{NA}).
#' @param first_mode Numeric value for the "First Mode" tPOD (optional; default is \code{NA}).
#' @param xrange A numeric vector of length 2 specifying the x-axis range (example is \code{c(0, 5)}) for zooming. default is NULL.
#'
#' @return A \code{ggplot2} object showing the BMD density distribution and optional vertical tPOD reference lines.
#'
#' @import ggplot2
#' @export

plot_BMD_tPOD_density <- function(BMD_values, lowest = NA,
                                         percentile = NA,
                                         mean_value = NA,
                                         accumulation = NA,
                                         first_mode = NA, xrange = NULL,
                                  title_label = "Density Distribution of BMD Values") {

  tpod_values <- c(
    "Lowest" = lowest,
    "Percentile" = percentile,
    "First Mode" = first_mode,
    "Accumulation" = accumulation,
    "Mean" = mean_value
  )

  # Define color mapping
  color_map <- c(
    "Lowest" = "#357D8A",
    "Percentile" = "orange",
    "First Mode" = "#235888",
    "Accumulation" = "#68A5D4",
    "Mean" = "#F7BD03"
  )

  idx = which(!is.na(tpod_values))

  if(length(idx)>0){

    tpod_values = tpod_values[idx]
    color_map = color_map[idx]

    line_data = data.frame(stat_name = names(tpod_values), stat_value = tpod_values, color = color_map)
    line_data$stat_name = factor(line_data$stat_name, levels = line_data$stat_name)

    # Create the plot
    p = ggplot(data.frame(BMD = BMD_values), aes(x = BMD)) +
      geom_density(fill = "gray", color = "darkgray", alpha = 0.6) +  # Density plot
      geom_vline(data = line_data, aes(xintercept = stat_value, color = stat_name),
                 linetype = "dashed", size = 1) +
      scale_color_manual(values = setNames(line_data$color, line_data$stat_name)) +
      theme_minimal(base_size = 15) +
      labs(title = title_label,
           x = "BMD",
           y = "Density",
           color = "tPOD") +  # Add a title for the legend
      theme(
        axis.title.x = ggplot2::element_text(size = 16),
        axis.title.y = ggplot2::element_text(size = 16),
        axis.text.x = ggplot2::element_text(size = 15),
        axis.text.y = ggplot2::element_text(size = 16),
        strip.text = ggplot2::element_text(size = 16),  # group_by label size
        legend.title = ggplot2::element_text(size = 16), # legend title
        legend.text = ggplot2::element_text(size = 14)
      )

  }else{
    p = ggplot(data.frame(BMD = BMD_values), aes(x = BMD)) +
      geom_density(fill = "gray", color = "darkgray", alpha = 0.6)+
      theme_minimal(base_size = 15) +
      labs(title = title_label,
           x = "BMD",
           y = "Density",
           color = "tPOD") +  # Add a title for the legend
      theme(
        axis.title.x = ggplot2::element_text(size = 16),
        axis.title.y = ggplot2::element_text(size = 16),
        axis.text.x = ggplot2::element_text(size = 15),
        axis.text.y = ggplot2::element_text(size = 16),
        strip.text = ggplot2::element_text(size = 16),  # group_by label size
        legend.title = ggplot2::element_text(size = 16), # legend title
        legend.text = ggplot2::element_text(size = 14)
      ) #+ coord_cartesian(xlim = c(0, 5))
  }

  if(!is.null(xrange)){
    p = p+ coord_cartesian(xlim = xrange)
  }
  return(p)
}




