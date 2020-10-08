#' Get labels from a dataset
#'
#' This function pulls out all of the labels from a dataset
#'
#' @param ds dataset
#' @return named vector with values of the column labels and names of column names
#' @examples
#' get_labels(mtcars)
#' @export

get_labels <- function(ds) {
  labs <- sapply(1:dim(ds)[2],FUN=function(x) attr(ds[[x]],"label"))
  names(labs) <- colnames(ds)

  labs
}

