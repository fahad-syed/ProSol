readTree <-
  function(file,restabs){
    library(ape)
    t <- read.tree(file)
    #Remove the "_species" names
    t$tip.label <- gsub("_species$","",t$tip.label)
    #Remove any node labels and brach lengths
    t$node.label <- NULL
    t$edge.length <- NULL
    common <- intersect(colnames(restabs$pfammat),t$tip.label)
    cat("readTree: From ",length(t$tip.label),"species in tree and",length(colnames(restabs$pfammat)),"species in results found",length(common),"common species.\n")
    #print(common)
    #Remove extra species from tree
    toremove <- setdiff(t$tip.label,common)
    t <- drop.tip(t,toremove)
    #plot(t)
    #Remove extra species from restabs
    toremove <- setdiff(colnames(restabs$pfammat),common)
    if (length(toremove) > 0) {
      restabs$pfam <- restabs$pfam[,-toremove]
      restabs$pfammat <- restabs$pfammat[,-toremove]
      restabs$clus <- restabs$clus[,-toremove]
      restabs$clusmat <- restabs$clusmat[,-toremove]
    }
    #Order restabs by tree
    restabs$pfam <- restabs$pfam[,c(colnames(restabs$pfam)[1:3],t$tip.label,grep("_IDs",colnames(restabs$pfam),value=T))]
    restabs$pfammat <- restabs$pfammat[,t$tip.label]
    restabs$clus <- restabs$clus[,c(colnames(restabs$clus)[1:2],t$tip.label,grep("_IDs",colnames(restabs$clus),value=T))]
    restabs$clusmat <- restabs$clusmat[,t$tip.label]
    out <- list(tree=t,restabs=restabs)
    invisible(out)
  }
