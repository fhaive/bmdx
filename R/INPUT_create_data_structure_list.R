#' this function convert the data in input into the format required for dose-dependent modelling
#' @param experimental_data a list of dataframe containing experimental data. Each row is a feature (e.g. gene) and each column is a sample
#' @param metadata a list of dataframe containing the metadata for the experimental data. Each row is a sample and the columns represent the different variables. A column for dose/concentration is required
#' @param sample_id_col a character specifying the name of the column containing the samples id
#' @param dose_id_col a character specifying the name of the column containing the doses/concentration
#' @param other_variables_id_col a vector of characters specifying the name of the column used to group the data
#' @param x a characters specifying the name of the x variable in the model. Default is dose.
#' @param y a characters specifying the name of the y variable in the model. Default is expr.
#' @return a dictionary containing the data frame for modelling. Dictionary Keys are n-uple specifying the experiment name, other variables of interests and feature names (e.g. drug, time, gene)
#' @export
create_data_structure = function(experimental_data,
                                      metadata,
                                      sample_id_col = "BARCODE",
                                      dose_id_col = "DOSE",
                                      other_variables_id_col = c("SACRI_PERIOD"),
                                      x = "dose",
                                      y = "expr"){

  experiment_names = names(experimental_data)

  for (mv in 1:length(metadata)) {
    rownames(metadata[[mv]]) = metadata[[mv]][,sample_id_col]
  }

  data_list = list()

  for (experiment in experiment_names) { # excel tabs
    print(experiment)
    combo_list = list()
    for (var in other_variables_id_col) { # per each variable
      combo_list[[var]] = unique(metadata[[experiment]][,var])
    }

    combinations = expand.grid(combo_list)

    pb = utils::txtProgressBar(min = 0, max = nrow(combinations), style = 3)
    for (row_idx in 1:nrow(combinations)) {

      is_sample_included = rep(TRUE, nrow(metadata[[experiment]]))
      for (col_idx in colnames(combinations)) {
        is_sample_included = is_sample_included & metadata[[experiment]][,col_idx] == combinations[row_idx,col_idx]
      }

      if (!any(is_sample_included)) { next() }

      samples_to_be_included = metadata[[experiment]][is_sample_included,sample_id_col]
      subset_experiment = as.matrix(experimental_data[[experiment]][,samples_to_be_included])
      doses_to_be_included = metadata[[experiment]][samples_to_be_included,dose_id_col]

      info = paste(as.vector(as.matrix(combinations[row_idx,])),collapse = "_")

      res = apply(subset_experiment,MARGIN = 1,FUN = function(riga){
        df = data.frame(doses_to_be_included,riga)
        colnames(df) = c(x,y)
        rownames(df) = NULL
        return(df)
      })

      names(res) = paste(experiment, info, names(res),sep = "_")

      data_list = append(data_list, res)
      utils::setTxtProgressBar(pb, row_idx)
    }
    close(pb)
  }
  return(data_list)
}
