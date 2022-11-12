#' Predict FEV1 based on Global Lungs Initiative
#'
#' This function computes predicted FEV1 from the ERS Global Lung Initiative function (Eur Resp J 2012 40:1324-43)
#'
#' @param age age in years (ages 3-95)
#' @param height height in <units>
#' @param race 1=African American (AA), North East Asian (NEA), South East Asian (SEA) or other (O)
#' @param sex 0=female/1=male
#' @param fev1 measured FEV1 (optional). If provided, Z-score and percent predicted will be computed
#' @return vector of predicted FEV1, coefficient of variation, skewness, lower limit normal,
#' Z-score (if fev1 is provided), percent predicted FEV1 (if fev1 is provided)
#' @keywords FEV1
#' @examples
#' @export


predict_FEV1 <- function(age, height, race, sex) {

  c


# SAS code example from Laura Colangelo
  # if sex=1 and race=5 then do;
  # Mspline = 0.3901 - 1.0579*(&age/100) + 1.4743*(&age/100)**2 - 2.1077*(&age/100)**3 - 0.1215*(&age/100)**4 +
  #   0.8873*(&age/100)**5;
  # M = exp(-10.3420 + 2.2196*log(&height) + 0.0574*log(&age) + Mspline);
  # Sspline = -1.6902 + 17.0986*(&age/100) - 68.1649*(&age/100)**2 + 127.1964*(&age/100)**3 - 109.6777*(&age/100)**4 +
  #   35.6832*(&age/100)**5;
  # S = exp(-2.3268 + 0.0798*log(&age) + Sspline);
  # L = 0.8866 + 0.0850*log(&age);
  # end;
  #
  # if sex=1 and race=4 then do;
  # Mspline = 0.3901 - 1.0579*(&age/100) + 1.4743*(&age/100)**2 - 2.1077*(&age/100)**3 - 0.1215*(&age/100)**4 +
  #   0.8873*(&age/100)**5;
  # M = exp(-10.3420 + 2.2196*log(&height) + 0.0574*log(&age) - 0.1589 + Mspline);
  # Sspline = -1.6902 + 17.0986*(&age/100) - 68.1649*(&age/100)**2 + 127.1964*(&age/100)**3 - 109.6777*(&age/100)**4 +
  #   35.6832*(&age/100)**5;
  # S = exp(-2.3268 + 0.0798*log(&age) + 0.1096 + Sspline);
  # L = 0.8866 + 0.0850*log(&age);
  # end;
  #
  # /*** Mspline for white and black women **/
  #   if sex=2 and race=5 then do;
  #   Mspline = 0.0552 + 1.6032*(&age/100) - 6.4855*(&age/100)**2 + 10.2741*(&age/100)**3 - 9.8646*(&age/100)**4 +
  #     3.8808*(&age/100)**5;
  #   M = exp(-9.6987 + 2.1211*log(&height) - 0.0270*log(&age) + Mspline);
  #   Sspline = -0.0822 + 1.4115*(&age/100) - 11.2797*(&age/100)**2 + 29.4613*(&age/100)**3 - 29.5597*(&age/100)**4 +
  #     10.4468*(&age/100)**5;
  #   S = exp(-2.3765 + 0.0972*log(&age) + Sspline);
  #   L = 1.1540 ;
  #   end;
  #
  #   if sex=2 and race=4 then do;
  #   Mspline = 0.0552 + 1.6032*(&age/100) - 6.4855*(&age/100)**2 + 10.2741*(&age/100)**3 - 9.8646*(&age/100)**4 +
  #     3.8808*(&age/100)**5;
  #
  #   M = exp(-9.6987 + 2.1211*log(&height) - 0.0270*log(&age) - 0.1484 + Mspline);
  #   Sspline = -0.0822 + 1.4115*(&age/100) - 11.2797*(&age/100)**2 + 29.4613*(&age/100)**3 - 29.5597*(&age/100)**4 +
  #     10.4468*(&age/100)**5;
  #   S = exp(-2.3765 + 0.0972*log(&age) + 0.1016 + Sspline);
  #   L = 1.1540 ;
  #   end;

}
