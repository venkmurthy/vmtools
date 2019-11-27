#' Calculate QUICKI
#'
#' This function calculates quantitative insulin sensitivity check index
#' Reference: Katz A, et al. Journal of Clinical Endocrinology & Metabolism 2000;85:2402-10
#' DOI: 10.1210/jcem.85.7.6661
#' @param insulin in μU/mL
#' @param glucose in mass units (mg/dL)
#' @return QUICKI as numeric value
#' @keywords QUICKI
#' @examples
#' calcQUICKI(15.0,100.0)
#'
#' glucose.vals <- seq(75,175,by=25)
#' insulin.vals <- seq(2,20, by=1)
#' test.dat <- expand.grid(insulin=insulin.vals,glucose=glucose.vals)
#'
#' test.dat$QUICKI <- calcQUICKI(insulin=test.dat$insulin,glucose=test.dat$glucose)
#' @export

calcQUICKI <- function(insulin,glucose) {
  return(1/(log(insulin) + log(glucose)))
}
