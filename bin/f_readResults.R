readResults <-
  function (fold="./"){
    rbypfam <- read.table(file=paste(fold,"InterproPfamOrganismStats.txt",sep='/'),stringsAsFactors=F,sep="\t",header=T)
    #Make to matrix for clustering
    names(rbypfam) <- gsub("\\.1","_IDs",names(rbypfam))
    idcols <- grep("_IDs",names(rbypfam))
    rbypfammat <- rbypfam[,-c(1:3,idcols)]
    rbypfammat<- as.matrix(rbypfammat)
    rownames(rbypfammat) <- paste(rbypfam[,'InterproID'],gsub("\\.\\d+$","",rbypfam[,'PfamId'],perl=T),sep='_')
    rownames(rbypfam) <- rownames(rbypfammat)
    #rbypfammat <- rbypfammat[,-1:2]
    
    statbypfam <- read.table(file=paste(fold,"InterproPfamMAD.txt",sep='/'),stringsAsFactors=F,sep="\t",header=T)
    rownames(statbypfam) <- paste(statbypfam[,'InterproID'],gsub("\\.\\d+$","",statbypfam[,'PfamId']),sep='_')
                                        #Get the cluster
    file <- list.files(path=fold,pattern="^PfamClusterOrganismSelectedOrgsStats_.*.txt",rec=T)
    remove <- grep("\\/ClusterPfamOrganismStats\\/",file)
    file <- file[-remove]
    #print(file)
    rbyclus <- read.table(paste(fold,file[1],sep='/'),header=1,stringsAsFactors=F,sep="\t")
    rbyclus <- rbyclus[-(1:nrow(rbyclus)),]
    for (i in 1:length(file)) {
      d <- read.table(paste(fold,file[i],sep='/'),header=1,stringsAsFactors=F,sep="\t")
      #Remove other than the searched for PFAM if they would happen to end up here.
      ipr <- gsub("/.*$","",file[i])
      #print(ipr)
      ipr <- grep(ipr,rownames(rbypfam),value=T)
      if (length(ipr) > 1) {
        stop(ipr,"has more than one matches\n")
      }
      pfam <- rbypfam[ipr,'PfamId']
      #print(dim(d))
      d <- d[d[,1] == pfam,]
      #     print(dim(d))
      #print(ipr)
      if (nrow(d) == 0) {
        warning(ipr," ",pfam," how no rows.\n")
      } else {
        d[,1] <- rep(ipr,nrow(d))
#        print(head(d))
      }
      rbyclus<- rbind(rbyclus,d)
    }

                                        #print(file)
    names(rbyclus) <- gsub("\\.1","_IDs",names(rbyclus))
    idcols <- grep("_IDs",names(rbyclus))
    rbyclusmat <- rbyclus[,-c(1:2,idcols)]
    rbyclusmat<- as.matrix(rbyclusmat)
    rownames(rbyclusmat) <- paste(gsub("\\.\\d+$","",rbyclus[,'PfamID']),rbyclus[,'ClusterId'],sep='_')
    rownames(rbyclus) <-  rownames(rbyclusmat)
    invisible(list(pfam=rbypfam,pfammat=rbypfammat,clus=rbyclus,clusmat=rbyclusmat,statpfam=statbypfam))
  }
