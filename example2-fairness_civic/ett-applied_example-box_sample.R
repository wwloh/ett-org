rm(list=ls())
# load the dataset
source("ett-applied_example-1prep.R")
# load the function to estimate the ETT
source("../help_funs-ett-dr_aipw.R") 

# step A1: fit propensity score model for IPW
fitPS <- glm(DISCRIMT1 ~ GENDERT1 + ETH_MCT1, 
             family = binomial("logit"), 
             data = da36561)
# step A2: obtain predicted propensity scores
PS.hat <- predict.glm(fitPS, type = "response")
# step A3: calculate IPW
da36561[, "IPW"] <- PS.hat / (1 - PS.hat)
da36561[da36561$DISCRIMT1==1L, "IPW"] <- 1.0




# step A4: outcome regression model for X==0 subgroup
fitY.X0 <- lm(FAIR_SB_T2 ~ GENDERT1 * ETH_MCT1, 
              data = da36561, 
              subset = DISCRIMT1==0)






# step A5:  calculate ETT estimator
ett.dr <- EstimateETT(
  data = da36561, # name of dataset
  trt.name = "DISCRIMT1", # name of discrimination variable (X)
  Y.name = "FAIR_SB_T2", # name of outcome variable (Y)
  regfit = fitY.X0 # fitted outcome model (among subset with X==0)
)
round(ett.dr$ett, 2)
# [1] -0.34

# step B2: outcome model with only main effect for treatment
fitY <- lm(FAIR_SB_T2 ~ DISCRIMT1 + GENDERT1 * ETH_MCT1, 
           data = da36561, 
           weights = IPW)
round(coef(fitY)["DISCRIMT1"], 2)
# DISCRIMT1 
# -0.34


# bootstrap
ONE.EST <- function(SEED) {
  check.ok <- FALSE
  while(!check.ok) {
    if (SEED==1L) {
      # observed data
      DATA <- da36561
    } else {
      # bootstrap sample with replacement
      DATA <- da36561[sample(nrow(da36561),nrow(da36561),replace=TRUE),]
    }
    # check for perfect separation
    fitPS <- glm(DISCRIMT1~GENDERT1+ETH_MCT1, family=binomial("logit"), data=DATA)
    check.ok <- all(abs(coef(fitPS))<10)
  }
  PS.hat <- predict.glm(fitPS, type="response")
  DATA[,"IPW"] <- PS.hat/(1-PS.hat)
  DATA[DATA$DISCRIMT1==1L,"IPW"] <- 1.0    
  fitY.X0 <- lm(FAIR_SB_T2~GENDERT1*ETH_MCT1, data=DATA, subset=DISCRIMT1==0)
  ett.dr <- EstimateETT(
    data=DATA, 
    trt.name="DISCRIMT1", 
    Y.name="FAIR_SB_T2", 
    regfit=fitY.X0)
  fitY <- lm(FAIR_SB_T2~DISCRIMT1+GENDERT1*ETH_MCT1, data=DATA, weights=IPW)
  ett.wls <- coef(fitY)["DISCRIMT1"]
  return(unlist(c("DR"=ett.dr,"WLS"=ett.wls)))
}



# ONE.EST: function that resamples observations with replacement,
# then calculates both ETT estimators for that resampled dataset
RES <- sapply(1:10000L, ONE.EST)
round(t(rbind("EST"=RES[,1], # point estimate using observed data
              apply(RES[,-1],1,quantile, probs=c(.025,.975)))),2)
#                 EST  2.5% 97.5%
# DR.ett        -0.34 -0.52 -0.16
# WLS.DISCRIMT1 -0.34 -0.52 -0.16
