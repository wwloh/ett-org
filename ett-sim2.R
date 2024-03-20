rm(list=ls())
library("data.table")
source("help_funs-ett-dr_aipw.R")

OneData <- function(N=100L, XL.relation="linear", YL.relation="linear") {
  L1 <- runif(N, min=-sqrt(3), max=sqrt(3))
  L2 <- rbinom(N, size=1, prob=0.1)
  # individual treatment effect
  U <- (-0.6) + runif(N, min=-0.3, max=0.3)
  # generate potential outcomes
  if (YL.relation=="linear") {
    Y0 <- L1 - L1*L2 + rnorm(N)
  } else {
    Y0 <- sqrt(abs(L1+0.5)) + rnorm(N)
  }
  Y1 <- Y0 + U
  
  if (XL.relation=="linear") {
    X.ast <- L1 + L1*L2 - 6*U
  } else {
    X.ast <- abs(L1)*(3-(L1^2)) + 2*(L2-2)*pmin(abs(L1)^(2*L1),2) - 6*U
  }
  PS <- exp(X.ast)/(1+exp(X.ast))
  ## check that positivity assumption holds
  check_positivity <- FALSE
  while(check_positivity==FALSE) {
    X <- rbinom(N, size=1, prob=PS)
    ## check that each stratum of L2 has both values of X to maintain positivity
    check_positivity <- all(table(L2,X)>0) && sum(X==0)>0 && sum(X==1)>0
  }
  Y <- (X==0)*Y0 + (X==1)*Y1
  
  # observed data
  DATA <- data.table(L1,L2,X,Y)
  setkey(DATA)
  
  # true effect values
  ATE <- mean(Y1-Y0)
  ETT <- mean(Y1[X==1] - Y0[X==1])
  
  rm(L1,L2,X,Y)
  
  # formulae for regression models
  form.Y <- as.formula(Y~X+L1*L2+L1^2)
  form.Y.X0 <- as.formula(Y~L1*L2+L1^2)
  form.PS <- as.formula(X~L1*L2+L1^2)
  
  # estimate ATE
  fit1 <- lm(form.Y, data=DATA)
  ATE.hat <- coef(fit1)["X"]
  rm(fit1)
  
  # estimate ETT
  
  # unweighted outcome model among untreated
  fitY.X0 <- lm(form.Y.X0, data=DATA, subset=X==0)
  ett.reg2 <- EstimateETT(
    data=DATA, 
    trt.name="X", 
    Y.name="Y", 
    regfit=fitY.X0, 
    outcome.only=TRUE,
    subgroup.effects=NULL)
  
  # IPW
  fitPS <- glm(form.PS, family=binomial("logit"), data=DATA)
  PS.hat <- predict.glm(fitPS, type="response")
  DATA[,"IPW"] <- PS.hat/(1-PS.hat)
  DATA[X==1L,"IPW"] <- 1.0
  
  # IPW estimator without outcome model
  ett.ipw1 <- EstimateETT(
    data=DATA, 
    trt.name="X", 
    Y.name="Y", 
    regfit=NULL, 
    outcome.only=FALSE,
    subgroup.effects=NULL)
  
  # Hajek estimator
  fit2 <- lm(Y~X, data=DATA, weights=IPW)
  ett.ipw2 <- coef(fit2)["X"]
  rm(fit2)
  
  # Morgan & Winship: outcome model with only main effect for treatment
  fitY <- lm(form.Y, data=DATA, weights=IPW)
  ett.dr1 <- coef(fitY)["X"]
  rm(fitY)
  
  # saturated outcome model among untreated with IPW
  ett.dr2 <- EstimateETT(
    data=DATA, 
    trt.name="X", 
    Y.name="Y", 
    regfit=fitY.X0, 
    outcome.only=FALSE,
    subgroup.effects=NULL)
  rm(fitY.X0)
  
  unlist(list(
    "ate"=ATE,
    "ate.hat"=ATE.hat,
    "bias.ate.hat"=(ATE.hat-ATE),
    "ett"=ETT,
    "ett.reg2"=ett.reg2$ett,
    "bias.ett.reg2"=(ett.reg2$ett-ETT),
    "ett.ipw1"=ett.ipw1$ett,
    "bias.ett.ipw1"=(ett.ipw1$ett-ETT),
    "ett.ipw2"=ett.ipw2,
    "bias.ett.ipw2"=(ett.ipw2-ETT),
    "ett.dr1"=ett.dr1,
    "bias.ett.dr1"=(ett.dr1-ETT),
    "ett.dr2"=ett.dr2$ett,
    "bias.ett.dr2"=(ett.dr2$ett-ETT)
  ))
}

SIMSETTINGS <- expand.grid(N=c(200L,5000L,80000L),
                           XL.relation=c("linear","nonlinear"),
                           YL.relation=c("linear","nonlinear"))

RES.LIST <- NULL
ptm <- proc.time()[3]
for (ss in 1:nrow(SIMSETTINGS)) {
  RES <- replicate(n=10000L,
                   OneData(N=SIMSETTINGS[ss,"N"], 
                           XL.relation=SIMSETTINGS[ss,"XL.relation"],
                           YL.relation=SIMSETTINGS[ss,"YL.relation"]))
  RES.LIST[[ss]] <- RES
  cat(ss,"setting done", round((proc.time()[3]-ptm)/60,1), "minutes so far\n")
}

save.image("ett-sim2.Rdata")

load("ett-sim2.Rdata")

BIAS <- do.call(rbind,lapply(RES.LIST, function(RES) 
  round(rowMeans(RES,na.rm=TRUE),2)
))

cbind(SIMSETTINGS,BIAS)

library("xtable")
print(xtable(cbind(SIMSETTINGS,BIAS),
             digits=c(0,0,0,0,rep(2,ncol(BIAS)))),
      include.rownames=FALSE)

q()