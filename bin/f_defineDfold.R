a <- commandArgs()
##print(a)
cat(paste("\nExecuting defineDfold.R\n\n"))
abbr <- NA
if(length(a) == 4) {
  ofold=a[length(a)]
} else {
  if(!exists("ofold"))  {
    ofold<-"."
  }
}
if(length(a) == 5) {
  ofold=a[length(a) -1 ]
  abbr =a[length(a)]
} else {
  if(!exists("ofold"))  {
    ofold<-"."
  }
}
pfold <- paste(ofold,"/PICS",sep="")
tfold <- paste(ofold,"/TABLES",sep="")
if (!is.na(abbr)) {
  pfold <- paste(pfold,"/",abbr,sep='')
  tfold <- paste(tfold,"/",abbr,sep='')
  if (!file.exists(pfold)) {
    dir.create(pfold, showWarnings = TRUE, recursive = FALSE, mode = "0777")
        cat("Created",pfold,"\n")
  }
  if (!file.exists(tfold)) {
    dir.create(tfold, showWarnings = TRUE, recursive = FALSE, mode = "0777")
        cat("Created",tfold,"\n")
  }
} else {
   if (!file.exists(pfold)) {
    dir.create(pfold, showWarnings = TRUE, recursive = FALSE, mode = "0777")
        cat("Created",pfold,"\n")
  }
  if (!file.exists(tfold)) {
    dir.create(tfold, showWarnings = TRUE, recursive = FALSE, mode = "0777")
        cat("Created",tfold,"\n")
  }
  
}


cat("Variables set by runall.pl:\n")
cat(paste("Output folder: ",ofold,"\n"))

#dfold <- "/mnt/msa1000-2/shared_data/Fungi/20111026EsaFungiModels"
#dfold <- "~/bigdata/20111102reco"
dfold <- gsub("/results.*","/data",ofold) 
cat(paste("Data folder: ",dfold,"\n"))
rfold <- paste(gsub("/results.*","/results",ofold),"/R-objects",sep='')
if (!file.exists(rfold)) {
    dir.create(rfold, showWarnings = TRUE, recursive = FALSE, mode = "0777")
        cat("Created",rfold,"\n")
  }
cat(paste("R-object folder:  ",rfold,"\n"))


cat(paste("Picture folder: ",pfold,"\n"))

cat(paste("Table folder: ",tfold,"\n"))

cat("Variables set by defineDfold.R:\n")
pfamfold <- "/u1/Blast/Databases/Pfam27.0"
cat(paste("PFAM folder: ",pfamfold,"\n"))
#pfamfold <- "/u1/Blast/Databases/Pfam27.0/"
#cat(paste("PFAM folder: ",pfamfold,"\n"))

#compRfold <- "~/proj/20120606reco/R-objects"
#compRfold <- "~/proj/20120724reco/R-objects"
#cat(paste("Background R-object folder (compRfold):  ",compRfold,"\n"))
#species <- 'Tree'
#cat(paste("Species is",species,"\n"))

