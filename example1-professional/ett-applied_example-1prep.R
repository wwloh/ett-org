rm(list=ls())
library("haven")
library("data.table")
da26782 <- read_sav(file="ICPSR_26782/DS0001/26782-0001-Data.sav")
da26782 <- data.table(da26782)
setkey(da26782)

# retain variables for analysis
to_keep <- c("GENDER","Q4ITNON","AGE","DEMO4","Q2",
             "Q204",
             "WRKSTRES"
             )
da26782 <- da26782[,..to_keep]
setkey(da26782)
rm(to_keep)

# exposure
da26782[, SUPFAIR.bin := as.integer(Q204)]
table(da26782[, list(SUPFAIR.bin,Q204)])
da26782[, SUPFAIR.bin := (SUPFAIR.bin>1)*1L]
da26782[, SUPFAIR := SUPFAIR.bin]
da26782[, SUPFAIR.bin := NULL]
da26782[, Q204 := NULL]

# outcome
da26782[, WRKSTRES := WRKSTRES/6]

# remove observations with missing covariates, exposure, or outcome
setkey(da26782)
da26782 <- da26782[!apply(
  apply(da26782[, list(GENDER,Q4ITNON,AGE,DEMO4,Q2,
                       SUPFAIR,WRKSTRES)],
        1,is.na),2,any)]
setkey(da26782)
nrow(da26782)
summary(da26782)

round(da26782[, mean(SUPFAIR)],3)

