makeSpeciesTable <-
  function(file,filter=NA){
    data <- read.delim(file,row.names=1,stringsAsFactors=F)[,1:7]
    colnames(data) <- gsub("\\."," ",colnames(data))
    ord <- order(data$Phylum,data$Subphylum,data$Class)
    data <- data[ord,]
    rownames(data) <- gsub(" ","_",data$Name)
    rownames(data) <- gsub("_species$","",rownames(data))
    if (length(filter) > 1) {
      cat("makeSpeciesTable: filtering for ",length(filter),"\n")
       filter <- gsub("_species$","",filter)
       data <- data[filter,]
    }
    invisible(data)
  }
