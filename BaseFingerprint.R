#libraries and entry requirements
library(prospectr)
library(rootSolve)
library(dplyr)
library(MonoInc)
library(DescTools)
library(caret)
library(readxl)
closest<-function(xv,sv){
  xv[which(abs(xv-sv)==min(abs(xv-sv)))] }
BD <- read_xlsx("Database1281.xlsx")
Thresholds <- read_xlsx("Thresholds.xlsx")
FP <- data.frame()

#Loop for all molecules 
 
for(i in 1:1281) {
  
  #Fp and csv file loading  
  d<- read.csv(file = sub("xyz", i, "xyz.csv"))
  colnames(d) <- c("V1","V2")
  
  if(d$V1[1] < 500){
    d <- d[min(which(d$V1 > 500)):length(d$V1),]
  }
  row.names(d) <- seq(1,length(d$V1))
  
  fk <-read.table(file = sub("xyz", i, "xyz.fp"), header = FALSE)
  
  #Bits thresholds 
  SG <- Thresholds
  SG$Bytes <- rep(0,101)
  for(ii in fk){
    SG$Bytes[ii] <- 1
  }
  
  
  #Space for modules of baseline and SG filter
  
  
  
  #Derivative of absorbance (used unless SG filter is used) 
  sg1 <- diff(d$V2,differences = 1)
  sg1 <- sg1 + rep(0.00000001, each = length(sg1))
  sg1 <- c(0.00000001,sg1)
  df1 <- d$V1
  
  #Zeros localization 
  mz <- c(uniroot.all(approxfun(df1,sg1), interval = c(min(df1), 700.99)),
          uniroot.all(approxfun(df1,sg1), interval = c(701, 1200.99)),
          uniroot.all(approxfun(df1,sg1), interval = c(1201, 1600.99)),
          uniroot.all(approxfun(df1,sg1), interval = c(1601, 2000.99)),
          uniroot.all(approxfun(df1,sg1), interval = c(2001, 2400.99)),
          uniroot.all(approxfun(df1,sg1), interval = c(2401, 2800.99)),
          uniroot.all(approxfun(df1,sg1), interval = c(2801, 3200.99)),
          uniroot.all(approxfun(df1,sg1), interval = c(3201, 3600.99)),
          uniroot.all(approxfun(df1,sg1), interval = c(3601, max(df1))))
  
  #Local minimums tracing
  if(which(df1==closest(df1,mz[1]))<= 2){
    mz <- mz[!mz %in% mz[1]]
  } 
  
  if(which(df1==closest(df1,mz[length(mz)]))>= as.integer(length(df1)-2)){
    mz <- mz[!mz %in% mz[length(mz)]]
  } 
  
  MZ <- c()
  
  for (nn in 1:length(mz)) {
    monopoint <- which(df1==closest(df1,mz[nn]))
    
    if(!monotonic(sg1[as.integer(monopoint-2):as.integer(monopoint+2)])){
      MZ <- c(MZ,mz[nn])
    }else{}
  }
  rm(mz,nn,monopoint)
  
  SG$Zeros <- NA
  
  for(rr in fk$V1){
    SG$Zeros[rr] <- sum(MZ > SG$Lower[rr] & MZ < SG$Upper[rr])  
  }
  
  #Fingerprint generation
  top <- c()
  
  for (rr in 1:101) {
    if(SG$Bytes[rr] == 0){
      top <- c(top, 0)
    }else if(SG$Bytes[rr] == 1 && SG$Zeros[rr] == 0){
      top <- c(top, 0)
    }else if(SG$Bytes[rr] == 1 && SG$Zeros[rr] > 0){
      top <- c(top, 1)
    }
  }
  #Aggregation of fingerprints
  top <- c(BD$logP[i],top)
  FP <- rbind(FP,top)
  
  rm(list=ls()[! ls() %in% c("closest", "FP","BD", "Thresholds")])  
}

#Column names and saving
fpp <- c()
for(i in 1:101){
  fpp <- c(fpp, sub("xx",i,"fpxx"))
}
colnames(FP) <- c("LogP", fpp)

write.csv(FP, file = "./FP/.csv", row.names = F)

 



#SG filter module
sg1 <- c(savitzkyGolay(X = c(rep(0,each = 50), d$V2[1:as.integer(max(which(d$V1<1500))+50)]), 
                       p = 3, w = 101, m = 1),
         savitzkyGolay(X = c(d$V2[as.integer(min(which(d$V1>1500))-250):length(d$V2)], rep(0,each = 250)), 
                       p = 3, w = 501, m = 1))
df1 <- d$V1

#Baseline module

sg1 <- rep(0, each = length(d$V2))

for(ii in 1:length(d$V2)) {
  if(d$V2[ii] <= 0.09)
    sg1[ii] <- 0
  else {
    sg1[ii] <- d$V2[ii] - 0.09
  }
}



