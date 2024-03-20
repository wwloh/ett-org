rm(list=ls())
library("data.table")
load("ICPSR_36561/DS0001/36561-0001-Data.rda")
da36561.0001 <- data.table(da36561.0001)
setkey(da36561.0001)
load("ICPSR_36561/DS0003/36561-0003-Data.rda")
da36561.0003 <- data.table(da36561.0003)
setkey(da36561.0003)

# retain only observations who filled in both surveys
da36561.0003$ID <- as.integer(as.character(da36561.0003$ID))
da36561 <- merge(da36561.0001,da36561.0003,by="ID")
rm(da36561.0001,da36561.0003)

# retain variables for analysis
to_keep <- c("GENDERT1","ETH_MCT1",
             "DISCRIMT1",
             "USFAIR_1T2","USFAIR_2T2","USFAIR_3T2"
             )
da36561 <- da36561[,..to_keep]
setkey(da36561)
rm(to_keep)

# exposure
## How often do you feel that you, personally, have been discriminated against for any reason?
da36561[, DISCRIMT1.bin := as.integer(DISCRIMT1)]
table(da36561[, list(DISCRIMT1.bin,DISCRIMT1)])
da36561[, DISCRIMT1.bin := (DISCRIMT1.bin>2)*1L]
da36561[, DISCRIMT1 := DISCRIMT1.bin]
da36561[, DISCRIMT1.bin := NULL]

# outcome
da36561[, USFAIR_1T2 := as.integer(USFAIR_1T2)]
da36561[, USFAIR_1T2 := as.integer(USFAIR_1T2)]
da36561[, USFAIR_1T2 := as.integer(USFAIR_1T2)]

setkey(da36561)
da36561[, FAIR_SB_T2 := mean(c(USFAIR_1T2,USFAIR_2T2,USFAIR_3T2), na.rm=TRUE),
        by=key(da36561)]
da36561[, c("USFAIR_1T2","USFAIR_2T2","USFAIR_3T2") := NULL]

# remove observations with missing covariates, exposure, or outcome
setkey(da36561)
da36561 <- da36561[!apply(
  apply(da36561[, list(GENDERT1,ETH_MCT1,
                       DISCRIMT1,FAIR_SB_T2)],
        1,is.na),2,any)]
setkey(da36561)
nrow(da36561)
summary(da36561)

round(da36561[, mean(DISCRIMT1)],3)
