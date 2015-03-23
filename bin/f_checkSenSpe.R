checkSenSpe <-
  function(path="/etc/ProSol/results/",start=12,ofold="./",print=T, plotAll=F){
    #TODO
    #Instead of plotting only those that have variation in sen or spe
    #plot the families that have variation in either so that
    #colors will match
    #consider adding legend
    #look at individual cases and try centralisation ->
    #would it fix min(abs(sen - spe))
    #plot these results

                                        #Read the data in
    file=grep("^IPR",list.files(path=path,pattern="^*PF.*PercentAverageSpecficitySensitivity",rec=T),value=T)
    td <- read.table(paste(path,file[1],sep=''),header=1,row.names=1)
    nOfIs <- nrow(td)
    spemat <- matrix(NA,nrow=0,ncol=nOfIs)
    senmat <- spemat
    for (f in 1:length(file)) {
      ipr <-  gsub("IPR0","",gsub("/.*$","",file[f]))
      d <- read.table(paste(path,file[f],sep=''),header=1,row.names=1)
      d<- d[order(rownames(d)),]
      spemat <- rbind(spemat,d[,2])
      colnames(spemat) <- rownames(d)
      senmat <- rbind(senmat,d[,1])
      colnames(senmat) <- rownames(d)
      rownames(spemat)[f] <- ipr
      rownames(senmat)[f] <- ipr
                                        #print(d)
                                        #print(ipr)
    }
                                        #    print(spemat)
                                        #Make output table
    otab <- cbind(senmat,spemat)
    colnames(otab)[1:nOfIs] <- paste("Sen",colnames(senmat),sep='')
    colnames(otab)[(nOfIs+1) : (2*nOfIs) ] <- paste("Spe",colnames(senmat),sep='')
    rownames(otab) <- gsub("^","IPR0",rownames(otab))
    otab <- cbind(otab,apply(senmat,1,var))
    colnames(otab)[ncol(otab)] <- 'Var of Sen'
    otab <- cbind(otab,apply(spemat,1,var))
    colnames(otab)[ncol(otab)] <- 'Var of Spe' 
    print(dim(otab))
                                        #Plot matrices
    #here you need to construct the filter first and combine it
    filt <- (!(apply(senmat,1,var) == 0)) | (!(apply(spemat,1,var) == 0))
    mats <- list('Sensitivity'=senmat,'Specificity'=spemat,'Sensitivity w/only variable'=senmat[filt,,drop=FALSE],'Specificity w/only variable'=spemat[filt,,drop=FALSE])
    if (print == T) {
      postscript(file=paste(ofold,"/",'SenSpeVsInf.eps',sep=''),hei=4*length(mats),wi=6,hor=F,paper='special')
    }
    layout(matrix(c(1:length(mats)), length(mats), 1, byrow = FALSE))
    
    for (i in names(mats)) {
                                        #Who has the maximum variation?
      mat <- mats[[i]]
      sds <- apply(mat,1,sd)
                                        #print(sds)
      names(sds) <- rownames(mat)
      msds <- sds[sds == max(sds)]
      cat(names(msds),"has the maximum variation for",i,"\n")  
      matplot(t(mat),type='l',ylab=i,xlab='10 x Inflation',axes=F)
      axis(1,at = 1:ncol(spemat), labels = colnames(spemat),las=2)
      axis(2)
                                        #print(senmat[names(sds),])
                                        #lines(y=mat[names(msds),],x=1:nOfIs,col='black')
      lines(y=apply(mat,2,mean),x=1:nOfIs,col='black',lwd=2)      
    }
    if (print == T) {
      dev.off()
    }
    #plot all individual pfams with variation
    noipr <- nrow(mats$'Sensitivity w/only variable')
    if (plotAll == T) {
      norows<- ceiling(noipr/3)
      postscript(file=paste(ofold,"/",'SenSpeVsInfAll.eps',sep=''),hei=norows*3,wi=6,hor=F,paper='special')
      lay <- matrix(c(1:(norows*3)), nrow=norows, ncol=3, byrow = TRUE)
      #print(lay)
      layout(lay)
      for (i in 1:noipr) {
        tit <- rownames(mats[['Sensitivity w/only variable']])[i]
        iprmat <- rbind(mats[['Sensitivity w/only variable']][i,],mats[['Specificity w/only variable']][i,])
        mi <- mean(iprmat)
        #get the point of min diff
        aipr <- abs(iprmat[1,] - iprmat[2,])
        #print(tit)
        #print(aipr)
        inf01 <- grep(min(aipr),aipr)
        inf01 <- inf01[1]
        #print(inf01)       
        #print(iprmat)
        #plot non-normalised
        matplot(t(iprmat),type='l',ylab="Sen & Spe",xlab='10 x Inflation',main=tit,col=c('grey', 'pink'),lty=1,axes=F)
        axis(1,at = 1:ncol(spemat), labels = colnames(spemat),las=2)
        axis(2)
        #print(rep(as.numeric(inf01),2))
        #print(c(mi*1.1,mi*0.9))
        #plot non normalised min diff
        lines(rep(as.numeric(inf01),2),c(mi*1.1,mi*0.9),col='green')
        #get the point of normalised min diff  
        iprmat <- t(scale(t(iprmat),scale=F))
        aipr <- abs(iprmat[1,] - iprmat[2,])
        inf02 <- grep(min(aipr),aipr)
        inf02 <- inf02[1]
        #print(aipr)
        #print(rep(as.numeric(inf02),2))
        #print(c(mi*1.1,mi*0.9))
        #plot normalised
        lines(rep(as.numeric(inf02),2),c(mi*1.1,mi*0.9),col='blue')
        #add average for plotting
        iprmat <- iprmat + mi
        #print(iprmat)
        lines(1:ncol(iprmat),iprmat[1,],col='black',lty=2)
        lines(1:ncol(iprmat),iprmat[2,],col='red',lty=2)

        
      }
      dev.off()
    }
    
                                        #What if I would use the largest I with identical to some early I rule
#    Is <- rep(NA,nrow(senmat))
#    names(Is) <- rownames(senmat)
#    for (i in rownames(senmat)) {
#      infl <- start
#      sen <- senmat[i,as.character(infl)]
#      ind <- grep(infl,colnames(senmat))
#      if (! sd(senmat[i,]) == 0) {
#        for (j in colnames(senmat)[ind:ncol(senmat)]) {          
#          jsen <- senmat[i,as.character(j)]
#          if (sen == jsen) {
#            infl <- j
#          } else {
#            break
#          }
#        }
#                                        #cat(sen,i,j,"Selected:",infl,"R:",senmat[i,],"\n")
#        Is[i] <- infl
#      }      
#    }
    #Use the normalised min diff rule:
     Is <- rep(NA,nrow(senmat))
    names(Is) <- rownames(senmat)
#    for (i in 1:nOfIs) {
        for (i in rownames(mats$'Specificity w/only variable')) {
         # rname <- rownames(senmat)[i]
#        tit <- rownames(mats[['Sensitivity w/only variable']])[i]
        iprmat <- rbind(mats[['Sensitivity w/only variable']][i,],mats[['Specificity w/only variable']][i,])
        iprmat <- t(scale(t(iprmat),scale=F))
        aipr <- abs(iprmat[1,] - iprmat[2,])
        inf02 <- grep(min(aipr),aipr)
        inf02 <- inf02[1]
      #cat(i,inf02,colnames(senmat)[inf02],"\n")
        Is[i] <- colnames(senmat)[inf02]
      }    
    otab <- cbind(otab,Is)
    colnames(otab)[ncol(otab)] <- 'I used'
    write.table(otab,file=paste(ofold,"/",'SenSpe.txt',sep=''),quote=F,sep='\t')
    out <- list(Is=Is,mats=mats,otab=otab)
    invisible(out)
  }
