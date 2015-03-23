`plot.PCA` <-
  function (pca,data,lxlim=c(0,0),lylim=c(0,0),xxlim=c(0,0),xylim=c(0,0),pcx=2,pcy=1,TOPp=4,TOPt=4,sbg='red',spch=21,TOPx=3,cex=1.2,htext=F,stext=NA,noside=T,cex.lab=1) {
    require(RColorBrewer)
     pch=21
    pch2=23
    
         layout(matrix(c(1,1,2,2,3,4,5,6), 4, 2, byrow = FALSE),widths=c(2,1))
    if (noside) {
      layout(matrix(1:2))
    }
    
    xxlim <- c(min(pca$x[,pcx]) + xxlim[1], max(pca$x[,pcx])+ xxlim[2] )
    xylim <- c(min(pca$x[,pcy]) + xylim[1], max(pca$x[,pcy])+ xylim[2] )
    plot(pca$x[,pcx],pca$x[,pcy],pch=spch,xlab=paste('Principal component',pcx,"with",summary(pca)$importance[2,pcx]*100,"% of variation",sep=' '),ylab=paste('Principal component',pcy, "with",summary(pca)$importance[2,pcy]*100,"% of variation",sep=' '),bg=sbg,col=sbg,bty='l',xlim=xxlim,ylim=xylim,cex=cex,cex.lab=cex.lab)
    rmax <- c(tail(names(sort(pca$x[,pcx])),1),head(names(sort(pca$x[,pcx])),1))
    rmax <- c(rmax,c(tail(names(sort(pca$x[,pcy])),1),head(names(sort(pca$x[,pcy])),1)))
    #cat(rmax,"\n")
    if (!is.na(stext[1]) ) {
      stext <- unique(c(rmax,stext))
      text(pca$x[stext,pcx],pca$x[stext,pcy],gsub("_"," ",stext),cex=cex)     
    } else {
      text(pca$x[rmax,pcx],pca$x[rmax,pcy],rmax,cex=cex)     
    }
 
    xmaxH <- head(names(sort(pca$rotation[,pcx])),TOPx)
     xmaxT <- tail(names(sort(pca$rotation[,pcx])),TOPx)
        ymaxH <- head(names(sort(pca$rotation[,pcy])),TOPx)
        ymaxT <- tail(names(sort(pca$rotation[,pcy])),TOPx)
    maxes <- c(xmaxH,xmaxT,ymaxH,ymaxT)
    #make index for maxes
    maxesi <- vector()
    for (i in maxes) {
      maxesi <- c(maxesi,grep(i,rownames(pca$x)))
    }
 #   cat("maxes",maxes,"\n",maxesi,"\n")
    lxlim <- c(min(pca$rotation[,pcx]) + lxlim[1], max(pca$rotation[,pcx])+ lxlim[2] )
    lylim <- c(min(pca$rotation[,pcy]) + lylim[1], max(pca$rotation[,pcy])+ lylim[2] )
    plot(pca$rotation[,pcx],pca$rotation[,pcy],xlim=lxlim,ylim=lylim,pch=pch,col='blue',bg='blue',bty='l',xlab=paste('Correlation with principal component',pcx,sep=' '),ylab=paste('Correlation with principal component',pcy,sep=' '),cex.lab=cex.lab)

    #this is now different from maxes and maxes would
    #probably be more informative
    pca.TOPp <- rev(names(tail(sort((pca$sdev[pcx]^2 * pca$rotation[,pcx]^2 + pca$sdev[pcy]^2 * pca$rotation[,pcy]^2 ) / (pca$sdev[pcx]^2 + pca$sdev[pcy]^2) ),TOPp)))
    pca.TOPt <- rev(names(tail(sort((pca$sdev[pcx]^2 * pca$rotation[,pcx]^2 + pca$sdev[pcy]^2 * pca$rotation[,pcy]^2 ) / (pca$sdev[pcx]^2 + pca$sdev[pcy]^2) ),TOPt)))
  #  cat(pca.TOPp,"\n")
    points(pca$rotation[pca.TOPp,c(pcx,pcy)],col='red',pch=pch,bg='red',cex=cex*1.7)
    text(pca$rotation[pca.TOPt,c(pcx,pcy)],labels=rownames(pca$rotation[pca.TOPt,]),pos=1,cex=cex)
# plot the same fox maxes
    points(pca$rotation[maxes,c(pcx,pcy)],col='orange',pch=pch,bg='orange')
    text(pca$rotation[maxes,c(pcx,pcy)],labels=rownames(pca$rotation[maxes,]),pos=1,cex=cex)
    
 #   hTOP <- head(pca.TOPp,2)
    hTOP <- rmax
    cat("length(hTOP)",length(hTOP),"\n",hTOP,"\n")
    hTIT <- c(paste('Most positive with PC',pcx,sep=''),paste('Most negative with PC',pcx,sep=''),paste('Most positive with PC',pcy,sep=''),paste('Most negative with PC',pcy,sep=''))
    
    #t1 <- data[rownames(pca$rotation),hTOP]
    t1 <- t(data[hTOP,rownames(pca$rotation)])
    cat("dim(t1)",dim(t1),"\n")
    if (!noside) {
      for (i in 4:1) {
        nona <- t1[,i]
        nona <- nona[!is.na(nona)]
       hxlim <- c(floor(min(nona) - median(nona)/2),ceiling(max(nona)+ median(nona)/2))
  #  cat('hxlim',hxlim,"\n")
#    htmp <- hist(nona,xlim=hxlim,main=paste(i,": ", names(t1)[i],sep=''),xlab=pca.TOPp[i])
        htmp <- hist(nona,xlim=hxlim,xlab=hTOP[i], main=hTIT[i],cex.lab=cex*1.2,ylab='')
#    yvec <- floor(seq(min(htmp$counts),max(htmp$counts),length.out=length(maxes)))
       yvec <- rep(mean(range(htmp$counts)),length(maxes))
  #  points(t1[,i],yvec,pch=spch,col=sbg,bg=sbg)
    #   cat(maxes,"\n")
        points(t1[maxes,hTOP[i]],yvec,pch=spch[maxesi],col=sbg[maxesi],bg=sbg[maxesi])
        if(htext==T) {
          text(t1[maxes,hTOP[i]],yvec,labels=maxes,srt=90,pos=3,cex=cex)
        }
     }
    }
    invisible(list('hTOP'=hTOP, 'pca.TOPp'= pca.TOPp, 'pca.TOPt'= pca.TOPt,'maxes'=maxes))
  }

