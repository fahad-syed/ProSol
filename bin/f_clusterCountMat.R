clusterCountMat <-
  function(data,nocl=c(11,30)){
   out <- list()
   for (i in nocl) {
 #         dataN <- t(scale(t(data)))
#     #remove NaN that come from dividing with 0
#     dataN[is.na(dataN)] <- 0
#      cl <- hclust(dist(dataN,method='euclidian'),method='complete')
#      cl$cluster <- cutree(cl,i)
#      out[[paste('hPearComp',i,sep='')]] <- cl
     cl <- kmeans(data, i, nstart =10,  algorithm = c("Hartigan-Wong"))
     out[[paste('kmeans',i,sep='')]] <- cl
   }
    invisible(out)
  }

