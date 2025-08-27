#' Read all sheets from an Excel file.
#'
#' This function reads all sheets from an Excel file specified by the \code{filename}.
#' It returns the data as a list of data frames, one for each sheet.
#'
#' @param filename The path to the Excel file.
#' @param tibble If \code{TRUE}, the output data frames will be converted to tibbles.
#'               Default is \code{FALSE}.
#' @param first_col_as_rownames If \code{TRUE}, the first column of each sheet will be used as row names.
#'                              Default is \code{FALSE}.
#' @param is_rnaseq_raw_count If \code{TRUE}, the data is assumed to be RNA-Seq raw counts, and it will
#'                            be converted to log2 counts using the \code{limma::voom} function.
#'                            Default is \code{FALSE}.
#'
#' @return A list of data frames, one for each sheet in the Excel file.
#' @export
#'

read_excel_allsheets <- function(filename, tibble = FALSE,
                                 first_col_as_rownames = FALSE,
                                 is_rnaseq_raw_count = FALSE,
                                 check_numeric = TRUE,
                                 na = "NA") {
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X){

    y = readxl::read_excel(filename, sheet = X,na = na)
    y = as.data.frame(y)
    if (first_col_as_rownames) {
      rownames(y) = y[,1]
      y = y[,-1]
    }

    if(check_numeric){ # forcing columns to be numerical
      for(i in 1:ncol(y)){
        y[,i] = as.numeric(y[,i])
      }
    }

    if (is_rnaseq_raw_count) {
      log2y = limma::voom(counts = as.matrix(y),design = NULL,normalize.method = "quantile")
      log2y = log2y$E
      rownames(log2y) = rownames(y)
      colnames(log2y) = colnames(y)
      y = log2y
    }
    y
  })

  if (!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  return(x)

}
