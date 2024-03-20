rm(list=ls())
# load the dataset
source("ett-applied_example-1prep.R")
# load the function to estimate the ETT
source("../help_funs-ett-dr_aipw.R") 

# step A1: fit propensity score model for IPW
fitPS <- glm(SUPFAIR ~ GENDER + Q4ITNON + AGE + DEMO4 + Q2,
             family = binomial("logit"), 
             data = da26782)
# step A2: obtain predicted propensity scores
PS.hat <- predict.glm(fitPS, type = "response")
# step A3: calculate IPW
da26782[, "IPW"] <- PS.hat / (1 - PS.hat)
da26782[da26782$SUPFAIR==1L, "IPW"] <- 1.0




# step A4: outcome regression model for X==0 subgroup
fitY.X0 <- lm(WRKSTRES ~ GENDER * Q4ITNON + AGE + I(AGE^2) + DEMO4 + Q2, 
              data = da26782, 
              subset = SUPFAIR==0)






# step A5:  calculate ETT estimator
ett.dr <- EstimateETT(
  data = da26782, # name of dataset
  trt.name = "SUPFAIR", # name of exposure variable (X)
  Y.name = "WRKSTRES", # name of outcome variable (Y)
  regfit = fitY.X0 # fitted outcome model (among subset with X==0)
)
round(ett.dr$ett, 2)
# [1] 0.88

# step B2: outcome model with only main effect for treatment
fitY <- lm(WRKSTRES ~ SUPFAIR + GENDER * Q4ITNON + AGE + I(AGE^2) + DEMO4 + Q2, 
           data = da26782, 
           weights = IPW)
round(coef(fitY)["SUPFAIR"], 2)
# SUPFAIR 
# 0.89


# bootstrap
ONE.EST <- function(SEED) {
  check.ok <- FALSE
  while(!check.ok) {
    if (SEED==1L) {
      # observed data
      DATA <- da26782
    } else {
      # bootstrap sample with replacement
      DATA <- da26782[sample(nrow(da26782),nrow(da26782),replace=TRUE),]
    }
    # check for perfect separation
    fitPS <- glm(SUPFAIR~GENDER+Q4ITNON+AGE+DEMO4+Q2, family=binomial("logit"), data=DATA)
    check.ok <- all(abs(coef(fitPS))<10)
  }
  PS.hat <- predict.glm(fitPS, type="response")
  DATA[,"IPW"] <- PS.hat/(1-PS.hat)
  DATA[DATA$SUPFAIR==1L,"IPW"] <- 1.0    
  fitY.X0 <- lm(WRKSTRES~GENDER*Q4ITNON+AGE+I(AGE^2)+DEMO4+Q2, data=DATA, subset=SUPFAIR==0)
  ett.dr <- EstimateETT(
    data=DATA, 
    trt.name="SUPFAIR", 
    Y.name="WRKSTRES", 
    regfit=fitY.X0)
  fitY <- lm(WRKSTRES~SUPFAIR+GENDER*Q4ITNON+AGE+I(AGE^2)+DEMO4+Q2, data=DATA, weights=IPW)
  ett.wls <- coef(fitY)["SUPFAIR"]
  return(unlist(c("DR"=ett.dr,"WLS"=ett.wls)))
}



# ONE.EST: function that resamples observations with replacement,
# then calculates both ETT estimators for that resampled dataset
RES <- sapply(1:10000L, ONE.EST)
round(t(rbind("EST"=RES[,1], # point estimate using observed data
              apply(RES[,-1],1,quantile, probs=c(.025,.975)))),2)
#                EST 2.5% 97.5%
# DR.ett      0.88 0.68  1.09
# WLS.SUPFAIR 0.89 0.69  1.09
