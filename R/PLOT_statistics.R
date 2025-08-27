
#' Plot Histogram and Density
#'
#' Generate a histogram and density plot visualization based on provided data and variables.
#'
#' @param mod_stats A data frame containing the data for plotting.
#' @param y_val The variable to be plotted on the y-axis (numeric).
#' @param color_by The categorical variable used for color-coding bars.
#' @param group_by The grouping variable for additional data segmentation.
#' @param group_by2 The second grouping variable for further segmentation.
#' @param filter_column The column name for filtering the data.
#' @param filter_by The values used for filtering the data.
#' @param alpha_th The transparency level of the density plot (ranges between 0-1).
#' @return A ggplot2 histogram and density plot visualization.
#' @export
plot_histogram = function(mod_stats, y_val = "BMD",
                          color_by = NULL,
                          group_by = NULL,
                          group_by2 = NULL,
                          filter_column = NULL,
                          filter_by = NULL,
                          alpha_th = 0.5){

  if (!all(c(y_val, group_by, group_by2, filter_column) %in% colnames(mod_stats))) {
    print("The selected columns are not present in the dataframe")
    return(make_empty_plot())
  }

  if (!is.null(color_by)) {
    if (class(mod_stats[,color_by]) != "factor") {
      print("The column selected as color has to be a factor!")
      return(make_empty_plot())
    }
  }

  if (class(mod_stats[,y_val]) != "numeric") {
    print("The selected column must contain numbers!")
    return(make_empty_plot())
  }

  if (alpha_th < 0 | alpha_th > 1) {
    print("The alpha value has ranges between 0-1")
    return(make_empty_plot())
  }

  DF = mod_stats

  if (!is.null(filter_column)) {
    DF = filter_df(DF, filter_column = filter_column,filter_by = filter_by)
    if (is.null(DF)) return(make_empty_plot())
  }

  p = ggplot2::ggplot(data = DF, ggplot2::aes_string(x = y_val, fill = color_by)) +
    ggplot2::geom_histogram(ggplot2::aes(y = ..density..)) +
    ggplot2::geom_density(alpha = alpha_th)


  if (!is.null(group_by)) {
    if (!is.null(group_by2)) {
      f = stats::formula(paste0(group_by2,"~",group_by,sep = ""))
    }else{
      f = stats::formula(paste0("~",group_by,sep = ""))
    }
    p = p + ggplot2::facet_grid(f)
  }

  return(p)
}



#' Plot Scatter Plot
#'
#' Generate a scatter plot visualization based on provided data and variables.
#'
#' @param mod_stats A data frame containing the data for plotting.
#' @param x_val The variable to be plotted on the x-axis.
#' @param y_val The variable to be plotted on the y-axis.
#' @param color_by The categorical variable used for color-coding points.
#' @param group_by The grouping variable for additional data segmentation.
#' @param group_by2 The second grouping variable for further segmentation.
#' @param filter_column The column name for filtering the data.
#' @param filter_by The values used for filtering the data.
#' @return A ggplot2 scatter plot visualization.
#' @export
plot_scatter = function(mod_stats, x_val = "BMDL", y_val = "BMD",
                        color_by = "Model",
                        group_by = NULL,
                        group_by2 = NULL,
                        filter_column = NULL,
                        filter_by = NULL){

  if (!all(c(x_val, y_val, color_by, group_by, group_by2, filter_column) %in% colnames(mod_stats))) {
    print("The selected columns are not present in the dataframe")
    return(make_empty_plot())

  }

  if (class(mod_stats[,x_val]) != "numeric" || class(mod_stats[,y_val]) != "numeric") {
    print("The selected column must contain numbers!")
    return(make_empty_plot())
  }

  if (class(mod_stats[,color_by]) != "factor") {
    print("The column selected as color has to be a factor!")
    return(make_empty_plot())
  }

  DF = mod_stats

  if (!is.null(filter_column)) {
    DF = filter_df(DF, filter_column = filter_column,filter_by = filter_by)
    if (is.null(DF)) return(make_empty_plot())
  }

  p = ggplot2::ggplot(DF, ggplot2::aes_string(x = x_val, y = y_val, color = color_by)) +
    ggplot2::geom_point(shape = 1)

  if (!is.null(group_by)) {
    if (!is.null(group_by2)) {
      f = stats::formula(paste0(group_by2,"~",group_by,sep = ""))
    }else{
      f = stats::formula(paste0("~",group_by,sep = ""))
    }
    p = p + ggplot2::facet_grid(f)
  }

  return(p)
}


#' this function allows to plot the histogram of one numeric column of the model statistics
#' @param mod_stats A data frame containing the data for plotting.
#' @param category The category variable to be used for the pie chart sectors.
#' @param group_by The grouping variable for additional data segmentation.
#' @param group_by2 The second grouping variable for further segmentation.
#' @param filter_column The column name for filtering the data.
#' @param filter_by The values used for filtering the data.
#' @return a ggplot object
#' @export
plot_pie_chart = function(mod_stats,
                          category = "Model",
                          group_by = NULL,
                          group_by2 = NULL,
                          filter_column = NULL,
                          filter_by = NULL){

  if (!all(c(group_by, group_by2, filter_column) %in% colnames(mod_stats))) {
    print("The selected columns are not present in the dataframe")
    return(make_empty_plot())
  }

  DF = mod_stats

  if (!is.null(filter_column)) {
    DF = filter_df(DF, filter_column = filter_column,filter_by = filter_by)
    if (is.null(DF)) return(make_empty_plot())
  }

  df = DF %>%
    dplyr::group_by_(.dots = c(category, group_by,group_by2)) %>%
    dplyr::summarise(Freq = n())

  df2 = df %>%
    dplyr::group_by_(.dots = c(group_by,group_by2)) %>%
    dplyr::mutate(Percentage = Freq/sum(Freq))

  p = ggplot2::ggplot(df2, ggplot2::aes(x = "",y=Percentage, fill = Model)) +
    ggplot2::geom_col() +
    ggplot2::coord_polar(theta = "y")


  if (!is.null(group_by)) {
    if (!is.null(group_by2)) {
      f = stats::formula(paste0(group_by2,"~",group_by,sep = ""))
    }else{
      f = stats::formula(paste0("~",group_by,sep = ""))
    }
    p = p + ggplot2::facet_grid(f)
  }


  return(p)
}


#' this function filters the dataframe of the models statistic
#' @param mod_stats the dataframe containing all the statistics
#' @param filter_column a vector of dataframe columns to be used for filtering
#' @param filter_by a list containing the selected values to be mantained after filtering.
#' Each position of the list correspond to one of the selected column and contains a vector of admissible values
#' #The length of filter_column and filter_by must be the same
#' @return the filtered dataframe
#' @export
filter_df = function(mod_stats,
                     filter_column,
                     filter_by){


  if (!is.null(filter_column)) { # if null, no filter is applied

    if (length(filter_column) != length(filter_by)) {
      print("The number of column to filter, and the list of filtering attributes have a different length")
      return(NULL)
    }

    if (!all(filter_column %in% colnames(mod_stats))) { #if any of the specified column is not present in the dataframe NULL is returned
      print("The selected column is not available")
      return(NULL)
    }

    for (i in 1:length(filter_column)) {
      fc = filter_column[i]
      idx = which(mod_stats[,fc] %in% filter_by[[i]])
      if (length(idx) == 0) {
        print("The selected values for the column are not present in the table")
        return(NULL)
      }
      mod_stats = mod_stats[idx,]
    }
  }

  mod_stats
}

#' Filter DataFrame by Column Values
#'
#' Filters a data frame based on specified columns and filter values.
#'
#' @param mod_stats A data frame containing model statistics.
#' @param filter_column A character vector specifying the columns to filter.
#' @param filter_by A list of character vectors containing filter values for each column.
#'
#' @return Returns a filtered data frame if filtering conditions are met; otherwise, NULL.
#' @export
filter_df_no = function(mod_stats,
                        filter_column,
                        filter_by){

  if (!is.null(filter_column)) { # if null, no filter is applied

    if (length(filter_column) != length(filter_by)) {
      print("The number of column to filter, and the list of filtering attributes have a different length")
      return(NULL)
    }

    if (!all(filter_column %in% colnames(mod_stats))) { #if any of the specified column is not present in the dataframe NULL is returned
      print("The selected column is not available")
      return(NULL)
    }

    for (i in 1:length(filter_column)) {
      fc = filter_column[i]
      idx = which(mod_stats[,fc] %in% filter_by[[i]])
      if (length(idx) == 0) {
        return(NULL)
      }
      mod_stats = mod_stats[idx,]
    }
  }

  mod_stats
}

#' Create Empty ggplot2 Plot
#' Creates an empty ggplot2 plot with a void theme and no x-label.
#' @return Returns an empty ggplot2 plot.
#' @export
make_empty_plot = function(){
  p = ggplot2::ggplot() +
    ggplot2::theme_void() +
    ggplot2::xlab(NULL)
  return(p)
}

#' Aggregate Rows by Time and Other Features
#'
#' Aggregates rows of a data frame based on time and other specified features.
#'
#' @param mod_stats A data frame containing model statistics.
#' @param gen_feat The general feature to aggregate.
#' @param first_feat The first feature used in aggregation.
#' @param group_by A character vector specifying the grouping features.
#' @param filter_column A character vector specifying the columns to filter.
#' @param filter_by A list of character vectors containing filter values for each column.
#'
#' @return Returns a ggplot2 plot of aggregated data or an empty plot if conditions are not met.
#' @export
aggregate_rows_time = function(mod_stats,
                               gen_feat = "Feature",
                               first_feat ="SACRI_PERIOD",
                               group_by = c("Experiment","DILI"),
                               filter_column = NULL,
                               filter_by = list(c("acetaminophen","bucetin"))){

  if (!all(c(group_by, filter_column) %in% colnames(mod_stats))) {
    print("The selected columns are not present in the dataframe")
    return(make_empty_plot())
  }

  if ( !is.null(filter_column)) {
    mod_stats = filter_df(mod_stats, filter_column = filter_column,filter_by = filter_by)
    if (is.null(mod_stats)) return(make_empty_plot())
  }

  L = list()
  for (i in c(first_feat,group_by)) {
    L[[i]] = mod_stats[,i]
  }
  df = stats::aggregate(mod_stats[,gen_feat], by = L, FUN = length)
  colnames(df) = c(first_feat,group_by,"counts")

  p <- ggplot2::ggplot(data = df, ggplot2::aes(x = factor(df[,first_feat],
                levels = unique(sort(as.numeric(as.vector(df[,first_feat])),decreasing = F))), y = df[,"counts"])) +
               ggplot2::geom_bar(stat = "identity")

  if (!is.null(group_by)) {
    if (length(group_by) == 1) {
      f = stats::formula(paste0("~",group_by, sep = ""))
    }
    if (length(group_by) == 2) {
      f = stats::formula(paste(group_by,collapse = "~"))
    }
    p = p + ggplot2::facet_grid(f)
  }

  p = p + ggplot2::scale_x_discrete(name = first_feat)
  p = p + ggplot2::scale_y_continuous(name = "Number of features")

  return(p)
}

make_empty_plot = function(){
  p = ggplot() +
    theme_void() +
    xlab(NULL)
  return(p)
}
