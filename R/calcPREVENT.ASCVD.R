#' Calculate PREVENT 10 year ASCVD risk scores
#'
#' This function calculates PREVENT 10 year ASCVD risk scores per recent AHA recommendations.
#' Reference: Khan SS, et al. Circulation 2023 online ahead of print
#' DOI: https://doi.org/10.1161/CIRCULATIONAHA.123.067626
#' @param age age in years
#' @param sex 0=female/1=male
#' @param tc total cholesterol in mg/dl
#' @param hdl HDL cholesterol in mg/dl
#' @param sbp Systolic blood pressure in mmHg
#' @param bptx treatment for blood pressure 0=FALSE/1=TRUE
#' @param smoking active smoking 0=FALSE/1=TRUE
#' @param dm diabetes mellitus 0=FALSE/1=TRUE
#' @param statin use of statin medication 0=FALSE/1=TRUE
#' @param egfr estimated glomerular filtration rate
#' @param female values which indicate female in the parameter sex, case-insensitive (default=0,"f","female")
#' @param male values which indicate male in the parameter sex, case-insensitive (default=1,"m","male")
#' @param bptx.true values which indicate BP med use, case-insensitive (default=1,"t","true","y","yes")
#' @param bptx.false values which indicate no BP med use, case-insensitive (default=0,"f","false","n","no")
#' @param smoking.true values which indicate active smoking, case-insensitive (default=1,"t","true","y","yes","active")
#' @param smoking.false values which indicate no active smoking, case-insensitive (default=0,"f","false","n","no","former","non-smoker","nonsmoker")
#' @param dm.true values which indicate diabetes, case-insensitive (default=1,"t","true","y","yes")
#' @param dm.false values which indicate no diabetes, case-insensitive (default=0,"f","false","n","no")
#' @param statin.true values which indicate statin use, case-insensitive (default=1,"t","true","y","yes")
#' @param statin.false values which indicate no statin use, case-insensitive (default=0,"f","false","n","no")
#' @return PREVENT 10 year ASCVD risk as numeric from 0 to 1
#' @keywords PREVENT
#' @examples
#' calcPREVENT.ASCVD(55,0,0,213,50,120,0,0,0)
#' calcPREVENT.ASCVD(55,1,0,213,50,120,0,0,0)
#' calcPREVENT.ASCVD(55,0,1,213,50,120,0,0,0)
#' calcPREVENT.ASCVD(55,1,1,213,50,120,0,0,0)
#'
#'
#' test.dat <- data.frame(age=rep(56:65,4),race=c(rep(0,10),rep(1,10),rep(0,10),rep(1,10)),
#'                        sex=c(rep(0,20),rep(1,20)),tc=rep(136:145,4),hdl=rep(46:55,4),sbp=rep(126:135,4),
#'                        bptx=rep(c(rep(0,5),rep(1,5)),4),smoking=rep(c(rep(0,5),rep("active",5)),4),
#'                        dm=rep(c(rep("No",5),rep(1,5)),4))
#'
#' test.dat$PREVENT.ASCVD <- calcPREVENT(test.dat$age,test.dat$race,test.dat$sex,test.dat$tc,test.dat$hdl,test.dat$sbp,
#'                         test.dat$bptx,test.dat$smoking,test.dat$dm)
#' @export

calcPREVENT.ASCVD <- function(age, sex, tc, hdl, sbp, bptx, smoking, dm, statin,
                        female=c(0,"f","female"),male=c(1,"m","male"),
                    bptx.true=c(1,"t","true","y","yes"), bptx.false=c(0,"f","false","n","no"),
                    smoking.true=c(1,"t","true","y","yes","active"),
                    smoking.false=c(0,"f","false","former","non-smoker","nonsmoker"),
                    dm.true=c(1,"t","true","y","yes"),dm.false=c(0,"f","false","n","no"),
                    statin.true=c(1,"t","true","y","yes"),statin.false=c(0,"f","false","n","no")) {


  # Initialize vectors
  lodds <- rep(NA,length(age))

  age2 <- age
  age2[age2>79] <- 79
  age2[age2<30] <- 30

  bprx <- rep(NA,length(bptx))
  bprx[sapply(bptx,tolower) %in% sapply(bptx.true,tolower)] <- 1
  bprx[sapply(bptx,tolower) %in% sapply(bptx.false,tolower)] <- 0
  smoke <- rep(NA,length(smoking))
  smoke[sapply(smoking,tolower) %in% sapply(smoking.true,tolower)] <- 1
  smoke[sapply(smoking,tolower) %in% sapply(smoking.false,tolower)] <- 0
  diab <- rep(NA,length(dm))
  diab[sapply(dm,tolower) %in% sapply(dm.true,tolower)] <- 1
  diab[sapply(dm,tolower) %in% sapply(dm.false,tolower)] <- 0
  lipidtx <- rep(NA,length(statin))
  lipidtx[sapply(statin,tolower) %in% sapply(statin.true,tolower)] <- 1
  lipidtx[sapply(statin,tolower) %in% sapply(statin.false,tolower)] <- 0

  # Women
  women <- sapply(sex,tolower) %in% sapply(female,tolower)

  lodds[women] <- -3.819975 + 0.719883 * (age2[women] - 55)/10 + 0.1176967 * (tc[women] - hdl[women] - 3.5) -
    0.151185 * (hdl[women] - 1.3)/0.3 - 0.0835358 * (pmin(sbp[women],110) - 110)/20 +
    0.3592852*(pmax(sbp[women], 110) - 130)/20 + 0.8348585*(diab[women]) + 0.4831078*(smoke[women]) +
    0.4864619*(pmin(egfr[women], 60) - 60)/-15 + 0.0397779*(pmax(egfr[women], 60) - 90)/-15 +
    0.2265309*(bprx[women]) - 0.0592374*(lipidtx[women]) - 0.0395762*(bprx[women])*(pmax(sbp[women],110) - 130)/20 +
    0.0844423*(lipidtx[women])*(tc[women] - hdl[women] - 3.5) -
    0.0567839*(age[women] - 55)/10 * (tc[women] - hdl[women] - 3.5) +
    0.0325692*(age[women] - 55)/10 * (hdl[women] - 1.3)/0.3 -
    0.1035985*(age[women] - 55)/10 * (pmax(sbp[women], 110) - 130) /20 -
    0.2417542*(age[women] - 55)/10 * (diab[women]) -
    0.0791142*(age - 55) /10 * (smoke[women]) - 0.1671492 * (age - 55) /10 * (pmin(eGFR, 60) - 60) / -15

  # Men
  men <- sapply(sex,tolower) %in% sapply(male,tolower)

  lodds[men] <- -3.500655 + 0.7099847*(age[men] - 55)/10 + 0.1658663*(tc[men] - hdl[men] - 3.5) -
    0.1144285*(hdl[men] - 1.3)/0.3 - 0.2837212*(pmin(sbp[men], 110) - 110)/20 +
    0.3239977*(pmax(sbp[men], 110) - 130)/20 + 0.7189597*(diab[men]) + 0.3956973*(smoke[men]) +
    0.3690075*(pmin(egfr[men], 60) - 60)/-15 + 0.0203619*(pmax(egfr[men], 60) - 90)/-15 +
    0.2036522*(bprx[men]) - 0.0865581*(lipidtx[men]) - 0.0322916*(bprx[men])*(pmax(sbp[men], 110) - 130)/20 +
    0.114563*(lipidtx[men])*(tc[men] - hdl[men] - 3.5) - 0.0300005*(age[men] - 55)/10 * (tc[men] - hdl[men] - 3.5) +
    0.0232747*(age[men] - 55)/10 * (hdl[men] - 1.3)/0.3 - 0.0927024*(age[men] - 55)/10 * (pmax(SBP, 110) - 130)/20 -
    0.2018525*(age[men] - 55)/10 * (diab[men]) - 0.0970527*(age[men] - 55)/10 * (smoke[men]) -
    0.1217081*(age[men] - 55)/10 * (pmin(egfr[men], 60) - 60)/-15

  prevent <- exp(lodds)/(1+exp(lodds))

  return(prevent)
}