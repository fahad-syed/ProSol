                                        #cat r_analyseClustering.R | R --slave --args /etc/ProSol/results/  
#source("f_defineDfold.R")
source("f_runFunction.R")
#library(sybil)
cat("\nExecuting r_analyseClustering.R\n\n")
a <- commandArgs()
fold <- a[4]
ofold <- paste(fold,"/PicsTables",sep="")
if (!file.exists(ofold)) {
  dir.create(ofold, showWarnings = TRUE, recursive = FALSE, mode = "0777")
  cat("Created",ofold,"\n")
}
#start = 10 x I-value to start the looking for a drop in sensitivity
start=12
if (length(a) == 5) {
  start <- a[5]
}
senspe <- paste(ofold,"/",'senspe.Rdata',sep='')
plotAll <-  TRUE
senspe <- runFunction(funcfile="f_checkSenSpe.R",outputfile=senspe,NamedListOfAddArgs=list(path='fold',start='start',ofold='ofold',plotAll='plotAll'))

#Load some data from a previous R-script in here 
#load(paste(rfold,"/03PCE.clustercols.Rdata",sep=''))
#Notice that plotAll plots also those which have no variation in either
#sensitivity or specificity and for which average inflation value
#was applied for.

