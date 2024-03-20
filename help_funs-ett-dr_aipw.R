EstimateETT <- function(data, # observed data
                        trt.name, # name of treatment variable
                        Y.name, # name of outcome variable
                        regfit=NULL, # fitted regression model(s)
                        outcome.only=FALSE, # use outcome model only
                        subgroup.effects=NULL # variables that may modify effect
                        ) {
  
  if (is.null(regfit)) {
    # must be IPW estimator if no regression model is specified
    outcome.only <- FALSE
  }
  D <- data.frame(data)
  TRT <- D[,trt.name]
  if (outcome.only==TRUE) {
    # subset with X==1
    D <- D[TRT==1,]
    TRT <- D[,trt.name]
    W <- 1 # for estimating the ETT
  } else {
    W <- (TRT-(1-TRT))*D[,"IPW"]
  }
  Y0.hat <- 0.0
  if (!is.null(regfit)) {
    # check if treatment is a predictor in fitted fitted regression model
    fitY.all <- trt.name %in% names(coef(regfit))
    if (fitY.all==TRUE) {
      D[,trt.name] <- 0
    }
    Y0.hat <- predict(regfit, newdata=D)
  }
  Ydiff <- D[, Y.name] - Y0.hat
  res <- NULL
  res[["ett"]] <- sum(W*Ydiff)/sum(TRT)
  if(!is.null(subgroup.effects)) {
    cett <- cett.names <- NULL
    subgroup.strata <- unique(D[,subgroup.effects,drop=FALSE])
    colnames(subgroup.strata) <- subgroup.effects
    row.names(subgroup.strata) <- NULL
    for (j in 1:nrow(subgroup.strata)) {
      in.subgroup <- apply(sapply(subgroup.effects, function(var.name) 
        D[,var.name] == subgroup.strata[j,var.name]),1,all)
      cett[[j]] <- sum((W*Ydiff)[in.subgroup])/sum(TRT[in.subgroup])
      cett.names[j] <- paste(sapply(1:length(subgroup.effects), function(i) 
        paste(c(subgroup.effects[i],subgroup.strata[j,i]),collapse="=")),
        collapse="__")
    }
    names(cett) <- cett.names
    res[["cett"]] <- cett
  }
  return(res)
}
