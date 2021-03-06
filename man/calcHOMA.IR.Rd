% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calcHOMA.IR.R
\name{calcHOMA.IR}
\alias{calcHOMA.IR}
\title{Calculate HOMA-IR}
\usage{
calcHOMA.IR(insulin, glucose)
}
\arguments{
\item{insulin}{in μU/mL}

\item{glucose}{in mass units (mg/dL)}
}
\value{
HOMA-IR as numeric value
}
\description{
This function calculates Homeostatic Model Assessment Insulin Resistance
Reference: Matthews DR, et al. Diabetologia 1995;28:412-9
DOI: 10.1007%2FBF00280883
}
\examples{
calcHOMA.IR(15.0,100.0)

glucose.vals <- seq(75,175,by=25)
insulin.vals <- seq(2,20, by=1)
test.dat <- expand.grid(insulin=insulin.vals,glucose=glucose.vals)

test.dat$HOMA.IR <- calcHOMA.IR(insulin=test.dat$insulin,glucose=test.dat$glucose)
}
\keyword{HOMA-IR}
