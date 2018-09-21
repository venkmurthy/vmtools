#' Calculate Pooled Cohorts Equations risk scores
#'
#' This function calculates Pooled Cohorts Equations risk scores per recent ACC/AHA recommendations.
#' Reference: Goff DC, et al. Circulation 2014;19:S49-S73
#' DOI: 10.1161/01.cir.0000437741.48606.98
#' @param age age in years
#' @param race white (non-black) vs. black
#' @param sex 0=female/1=male
#' @param tc total cholesterol in mg/dl
#' @param hdl HDL cholesterol in mg/dl
#' @param sbp Systolic blood pressure in mmHg
#' @param bptx treatment for blood pressure 0=FALSE/1=TRUE
#' @param smoking active smoking 0=FALSE/1=TRUE
#' @param dm diabetes mellitus 0=FALSE/1=TRUE
#' @param white value which indicates white/non-black the parameter race (default=0)
#' @param black value which indicates blacks in the parameter race (default=1)
#' @param female value which indicates female in the parameter sex (default=0)
#' @param male value which indicates male in the parameter sex (default=1)
#' @return PCE risk as numeric from 0 to 1
#' @keywords PCE
#' @examples
#' calcPCE(55,0,0,213,50,120,0,0,0)
#' calcPCE(55,1,0,213,50,120,0,0,0)
#' calcPCE(55,0,1,213,50,120,0,0,0)
#' calcPCE(55,1,1,213,50,120,0,0,0)
#'
#'
#' test.dat <- data.frame(age=rep(56:65,4),race=c(rep(0,10),rep(1,10),rep(0,10),rep(1,10)),
#'                        sex=c(rep(0,20),rep(1,20)),tc=rep(136:145,4),hdl=rep(46:55,4),sbp=rep(126:135,4),
#'                        bptx=rep(c(rep(0,5),rep(1,5)),4),smoking=rep(c(rep(0,5),rep(1,5)),4),
#'                        dm=rep(c(rep(0,5),rep(1,5)),4))
#'
#' test.dat$PCE <- calcPCE(test.dat$age,test.dat$race,test.dat$sex,test.dat$tc,test.dat$hdl,test.dat$sbp,
#'                         test.dat$bptx,test.dat$smoking,test.dat$dm)
#' @export

calcPCE <- function(age, race, sex, tc, hdl, sbp, bptx, smoking, dm, white=0, black=1, female=0, male=1) {

  # Initialize vector
  pce <- rep(NA,length(age))

  # White women
  ww <- sex==female & race==white
  pce[ww] <- 29.18 +
    (-29.799 * log(age[ww])) + (4.884 * log(age[ww])^2) +
    (13.540 * log(tc[ww])) + (-3.114 * log(age[ww]) * log(tc[ww])) + (-13.578 * log(hdl[ww])) +
    (3.149 * log(age[ww]) * log(hdl[ww])) + (2.019 * log(sbp[ww]) * bptx[ww]) +
    (1.957 * log(sbp[ww]) * as.numeric(bptx[ww]==0)) + (7.574 * smoking[ww]) +
    (-1.665 * log(age[ww]) * smoking[ww]) + (0.661 * dm[ww])

  # Black women
  bw <- sex==female & race==white
  pce[bw] <- -86.61 +
    (17.114 * log(age[bw])) + (0.940 * log(tc[bw])) + (-18.920 * log(hdl[bw])) +
    (4.475 * log(age[bw]) * log(hdl[bw])) + (29.291 * log(sbp[bw]) * bptx[bw]) +
    (-6.432 * log(age[bw]) * log(sbp[bw]) * bptx[bw]) +
    (27.820 * log(sbp[bw]) * as.numeric(bptx[bw]==0)) +
    (-6.087 * log(age[bw]) * log(sbp[bw]) * as.numeric(bptx[bw]==0)) +
    (0.691 * smoking[bw]) + (0.874 * dm[bw])

  # White men
  wm <- sex==male & race==white
  pce[wm] <- -61.18 +
    (12.344 * log(age[wm])) + (11.853 * log(tc[wm])) + (-2.664 * log(age[wm]) * log(tc[wm])) +
    (-7.990 * log(hdl[wm])) + (1.769 * log(age[wm]) * log(hdl[wm])) + (1.797 * log(sbp[wm]) * bptx[wm]) +
    (1.764 * log(sbp[wm]) * as.numeric(bptx[wm]==0)) + (7.837 * smoking[wm]) +
    (-1.795 * log(age[wm]) * smoking[wm]) + (0.658 * dm[wm])

  # Black men
  bm <- sex==white & race==black
  pce[bm] <- -19.54 +
    (2.469 * log(age[bm])) + (0.302 * log(tc[bm])) +
    (-0.307 * log(hdl[bm])) + (1.916 * log(sbp[bm]) * bptx[bm]) +
    (1.809 * log(sbp[bm]) * as.numeric(bptx[bm]==0)) + (0.549 * smoking[bm]) +
    (0.645 * dm[bm])

  pce <- exp(pce)
  pce[ww] <- 1 - 0.9665 ^ pce[ww]
  pce[bw] <- 1 - 0.9533 ^ pce[bw]
  pce[wm] <- 1 - 0.9144 ^ pce[wm]
  pce[bm] <- 1 - 0.8954 ^ pce[bm]

  pce
}

