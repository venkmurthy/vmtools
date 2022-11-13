#' Predict FEV1 based on Global Lungs Initiative
#'
#' This function computes predicted FEV1 from the ERS Global Lung Initiative function (Eur Resp J 2012 40:1324-43)
#'
#' @param age age in years (ages 25-95)
#' @param height height in centimeters
#' @param race 1=African American (AA), 2=North East Asian (NEA), 3=South East Asian (SEA) or 4=other (O)
#' @param sex 0=female/1=male
#' @param fev1 measured FEV1 (optional). If provided, Z-score and percent predicted will be computed
#' @return vector of predicted FEV1, coefficient of variation, skewness, lower limit normal,
#' Z-score (if fev1 is provided), percent predicted FEV1 (if fev1 is provided)
#' @keywords FEV1
#' @examples
#'
#' age <- c(4.8, 12.2, 53, 39.1)
#' height <- c(107, 152, 175, 165)
#' race <- c("o", "aa", "SEA", "SEA")
#' sex <- c("M", "M", "M", "F")
#' fev1 <- c(0.8, 2.405, 2.410, 2.210)
#'
#' predict_FEV1(age, height, race, sex, fev1)
#'
#' @export


predict_FEV1 <- function(age, height, race, sex, fev1=NA) {

  # Make sure length of all variables are equal
  stopifnot(all(sapply(list(age, height, race), function(x) length(x) = length(sex))))
  stopifnot(is.na(fev1) | (length(fev1)==length(age)))

  # Check if age is in the right range
  age <- ifelse(age <25 | age >95, NA, age)
  female <- ifelse(tolower(as.character(sex)) %in% c("0","f","female"),TRUE,FALSE)

  # Set up matrix of x-variables
  x_mat <- data.frame(intercept=1,
                      ln_height=log(height),
                      ln_age=log(age),
                      aa=ifelse(tolower(as.character(race)) %in% c("1","aa","african-american"),1,0),
                      nea=ifelse(tolower(as.character(race)) %in% c("2","nea","N east asian"),1,0),
                      sea=ifelse(tolower(as.character(race)) %in% c("3","sea","s east asian"),1,0),
                      o=ifelse(tolower(as.character(race)) %in% c("4","o","other"),1,0),
                      age.ge.25=ifelse(age>=25,1,0),
                      spline.age1=ifelse(age>=25,age/100,0),
                      spline.age2=ifelse(age>=25,(age/100)^2,0),
                      spline.age3=ifelse(age>=25,(age/100)^3,0),
                      spline.age4=ifelse(age>=25,(age/100)^4,0),
                      spline.age5=ifelse(age>=25,(age/100)^5,0)) |> as.matrix()

  # Set up matrix for M coefficients
  m_coef_mat <- data.frame(intercept=ifelse(female,-9.6987,-10.3420),
                           ln_height=ifelse(female,2.1211,2.2196),
                           ln_age=ifelse(female,-0.0270,0.0574),
                           aa=ifelse(female,-0.1484,-0.1589),
                           nea=ifelse(female,-0.0149,-0.0351),
                           sea=ifelse(female,-0.1206,-0.0881),
                           o=ifelse(female,0,0),
                           age.ge.25=ifelse(female,0.0552,0.3901),
                           spline.age1=ifelse(female,1.6032,-1.0579),
                           spline.age2=ifelse(female,-6.4855,1.4743),
                           spline.age3=ifelse(female,10.2741,-2.1077),
                           spline.age4=ifelse(female,-9.8646,-0.1215),
                           spline.age5=ifelse(female,3.8808,0.8873)) |> as.matrix()

  # Set up matrix for S coefficients
  s_coef_mat <- data.frame(intercept=ifelse(female,-2.3765,-2.3268),
                           ln_height=ifelse(female,0,0),
                           ln_age=ifelse(female,0.0972,0.0798),
                           aa=ifelse(female,0.1016,0.1096),
                           nea=ifelse(female,-0.0109,-0.3973),
                           sea=ifelse(female,0.0733,0.0327),
                           o=ifelse(female,0,0),
                           age.ge.25=ifelse(female,-0.0825,-1.6902),
                           spline.age1=ifelse(female,1.4104,17.0986),
                           spline.age2=ifelse(female,-11.2699,-68.1649),
                           spline.age3=ifelse(female,29.4400,127.1964),
                           spline.age4=ifelse(female,-29.5505,-109.6777),
                           spline.age5=ifelse(female,10.4405,35.6832)) |> as.matrix()

  # Set up matrix for L coefficients
  l_coef_mat <- data.frame(intercept=ifelse(female,1.1540,0.8866),
                           ln_height=ifelse(female,0,0),
                           ln_age=ifelse(female,0,0.0850),
                           aa=ifelse(female,0,0),
                           nea=ifelse(female,0,0),
                           sea=ifelse(female,0,0),
                           o=ifelse(female,0,0),
                           age.ge.25=ifelse(female,0,0),
                           spline.age1=ifelse(female,0,0),
                           spline.age2=ifelse(female,0,0),
                           spline.age3=ifelse(female,0,0),
                           spline.age4=ifelse(female,0,0),
                           spline.age5=ifelse(female,0,0))

  # Compute M, S, L, LLN, percent predicted, Z-score
  M <- exp(rowSums(x_mat * m_coef_mat))
  S <- exp(rowSums(x_mat * s_coef_mat))
  L <- rowSums(x_mat * l_coef_mat)
  LLN <- exp(log(1-1.644*L*S)/L + log(M))
  pct_pred <- (fev1/M) * 100
  Z <-((fev1/M)^L -1)/(L*S)

  # Data frame for return
  data.frame(fev1_pred=M, cv=S, skew=L, fev1_lln=LLN, fev1_pct_pred=pct_pred, fev1_Z=Z)
}

age <- c(4.8, 12.2, 53, 39.1)
height <- c(107, 152, 175, 165)
race <- c("o", "aa", "SEA", "SEA")
sex <- c("M", "M", "M", "F")
fev1 <- c(0.8, 2.405, 2.410, 2.210)

predict_FEV1(age, height, race, sex, fev1)
