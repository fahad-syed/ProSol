getClusterNosFromPCA <-
  function(pca,proplimits=c(0.9,0.95,0.99)){
    clusterNos <- vector()
    for (i in proplimits) {
      cumsum<- summary(pca)$importance[3,]
      PC <- names(cumsum)[cumsum > i]
      PC <- PC[1]
      PC <- as.numeric(gsub("PC","",PC))
      clusterNos <- c(clusterNos,PC)
      cat("getClusterNosFromPCA:",i,PC,"\n")
    }
    names(clusterNos) <- proplimits
    invisible(clusterNos)
  }
