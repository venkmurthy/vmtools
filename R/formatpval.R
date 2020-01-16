#' Format p-value for pretty printing
#'
#' @param p vector of p-values
#' @param include.p if TRUE then the characters "p=" or "p<" are prepended (default:TRUE)
#' @param na.string characters to output if p-value is NA (default: "NA")
#' @examples
#' p <- c(0.5,0.499999,0.004,0.0003,1e-10,NA)
#' format.pval(p)
#' @export

format.pval <- function(p,include.p=TRUE,na.string="NA") {
  if (include.p==TRUE) {
    ifelse(is.na(p),na.string,
           ifelse(p<0.0001,"p<0.0001",
                  ifelse(p<0.001,sprintf("p=%0.4f",p),
                         ifelse(p<0.01,sprintf("p=%0.3f",p),sprintf("p=%0.2f",p)))))
  } else {
    ifelse(is.na(p),na.string,
           ifelse(p<0.0001,"<0.0001",
                  ifelse(p<0.001,sprintf("%0.4f",p),
                         ifelse(p<0.01,sprintf("%0.3f",p),sprintf("%0.2f",p)))))
  }
}
