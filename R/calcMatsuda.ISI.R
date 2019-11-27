#' Calculate Matsuda ISI
#'
#' This function calculates Matsuda insulin sensitivity index from 2 hour oral glucose tolerance test
#' Reference: Matsuda M and DeFronzo RA. Diabetes Care 1999; 22:1462-70
#' DOI: 10.2337/diacare.22.9.1462
#' @param fasting.insulin in μU/mL
#' @param fasting.glucose in mass units (mg/dL)
#' @param mean.insulin in μU/mL
#' @param mean.glucose in mass units (mg/dL)
#' @return Matsuda ISI as numeric value
#' @keywords Matsuda ISI
#' @examples
#' calcMatsuda.ISI(15.0,100.0,30,150)
#'
#' fasting.glucose.vals <- seq(75,175,by=25)
#' fasting.insulin.vals <- seq(2,20, by=1)
#' mean.glucose.vals <- seq(75,175,by=25)
#' mean.insulin.vals <- seq(2,20, by=1)
#' test.dat <- expand.grid(fasting.insulin=insulin.vals,fasting.glucose=glucose.vals,
#'                         mean.insulin=mean.insulin.vals,mean.glucose=mean.glucose.vals)
#'
#' test.dat$Matsuda.ISI <- calcMatsuda.ISI(fasting.insulin=test.dat$fasting.insulin,
#'                                         fasting.glucose=test.dat$fastingglucose,
#'                                         mean.insulin=test.dat$mean.insulin,
#'                                         mean.glucose=test.dat$mean.glucose)
#' @export

calcMatsuda.ISI <- function(fasting.insulin,fasting.glucose,mean.insulin,mean.glucose) {
  return(10000/sqrt(fasting.insulin*fasting.glucose*mean.insulin*mean.glucose))
}
