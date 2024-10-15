#' Fast read wrapper for multiple objects using qs
#'
#' This reads output of mqsave which wraps objects in a list to allow saving multiple objects using qs
#'
#' @param objects vector containing names of objects to save
#' @param fname filename
#' @param ... additional parameters passed to qsave
#' @return returns result of qsave call from qs module
#' @keywords serialization fast save
#' @examples
#' calcHOMA.IR(15.0,100.0)
#'
#' test.data <- data.frame(var1=rnorm(1000), var2=rnorm(1000))
#' mqsave(objects="test.data",fname="test.qs")
#' rm(test.data)
#' mqread(fname="test.qs")
#' @export

mqread <- function(fname, ...) {
  # Call qread
  object.list <- qs::qread(file=fname)

  # Unlist objects
  for(x in names(object.list)) assign(x,object.list[[x]],envir=parent.frame())
}
