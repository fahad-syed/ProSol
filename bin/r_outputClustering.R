                                        #cat r_outputClustering.R | R --slave --args /etc/ProSol/results/2014-02-28_for_all_sequences
#source("f_defineDfold.R")
source("f_runFunction.R")
#library(sybil)
cat("\nr_outputClustering.R\n\n")
a <- commandArgs()
fold <- a[4]
ofold <- paste(fold,"/PicsTables",sep="")
if (!file.exists(ofold)) {
  dir.create(ofold, showWarnings = TRUE, recursive = FALSE, mode = "0777")
  cat("Created",ofold,"\n")
}

#Read in InterproPfamMAD.txt and InterproPfamOrganismStats.txt
#Read in the PfamClusterOrganismSelectedOrgsStats_??.txt and concatenate
restabs <- paste(ofold,"/",'restabs.Rdata',sep='')
restabs <- runFunction(funcfile="f_readResults.R",outputfile=restabs,NamedListOfAddArgs=list(fold='fold'))


#Read in tree
#take the intersect of table and tree
#Order results columns by tree
treerestabs <- paste(ofold,"/",'treerestabs.Rdata',sep='')
treefile <- "../data/originaltree.nwk"
treerestabs <- runFunction(funcfile="f_readTree.R",outputfile=treerestabs,NamedListOfAddArgs=list(file='treefile',restabs='restabs'))
tree <- treerestabs$tree
restabs <- treerestabs$restabs

#Read in fungalspecies.txt
species <- paste(ofold,"/",'species.Rdata',sep='')
speciesfile <- "../data/SpeciesTable.txt"
species <- runFunction(funcfile="f_makeSpeciesTable.R",outputfile=species,NamedListOfAddArgs=list(file='speciesfile',filter='tree$tip.label'))

#PCA
source("f_plot.PCA.R")
#Highlight species of interest
stext=c("Calcarisporium","Scopulariopsis","Pestalotiopsis")
#Leave Magnaporthe_grisea  and IPR001138 and IPR007219 
data <- t(restabs$pfammat[-which(rownames(restabs$pfammat) %in% c("IPR001138_PF00172","IPR007219_PF04082")),-which(colnames(restabs$pfammat) %in% c('Magnaporthe_grisea','Rhizopus_oryzae'))])
pca <- prcomp(data,retx=T)
postscript(file=paste(ofold,"/pcapfam1.eps",sep=''),horizontal=F,paper='special',width=7,hei=10)
data.plot.PCA <- plot.PCA(pca=pca,data=data,sbg=species[rownames(data),'Color'],spch= 21,lylim=c(-0.05,0.05),lxlim=c(-0.05,0.05),xxlim=c(-10,10),xylim=c(0,0),cex=0.6,cex.lab=0.8,stext=stext)
dev.off()

#Leave Magnaporthe_grisea &  Rhizopus_oryzae out from this data set as it is too different
data <- t(restabs$pfammat[,-which(colnames(restabs$pfammat) %in% c('Magnaporthe_grisea','Rhizopus_oryzae'))])

pca <- prcomp(data,retx=T)
postscript(file=paste(ofold,"/pcapfam2.eps",sep=''),horizontal=F,paper='special',width=7,hei=10)
data.plot.PCA <- plot.PCA(pca=pca,data=data,sbg=species[rownames(data),'Color'],spch= 21,lylim=c(-0.05,0.05),lxlim=c(-0.05,0.05),xxlim=c(-10,10),xylim=c(0,0),cex=0.6,cex.lab=0.8,stext=stext)
dev.off()
data <- t(data)

#Cluster
#Get number of clusters
source("f_getClusterNosFromPCA.R")
proplimits <- c(0.999,0.9995)
nocl<- getClusterNosFromPCA(pca,proplimits=proplimits)
clusters <- paste(ofold,"/",'clusters.Rdata',sep='')
clusters <- runFunction(funcfile="f_clusterCountMat.R",outputfile=clusters,NamedListOfAddArgs=list(data='data',nocl='nocl'))
#Cut matrix -> this does not work
#dpfammat <- paste(ofold,"/",'dpfammat.Rdata',sep='')
#breaks <- c(0,1:5,max(restabs$pfammat))
#labels <- c(0:5)
#dpfammat <- runFunction(funcfile="f_discretisemat.R",outputfile=dpfammat,NamedListOfAddArgs=list(data='restabs$pfammat',breaks='breaks',labels='labels'))
#data <- dpfammat$data
data <- restabs$pfammat
data[data > 5 ] <- 5

#Get interpro names
#grep "interpro id" /mnt/msa1000-2/interpro/iprscan/data/interpro.xml |  perl -pe "s/^.*id=//" | perl -pe "s/pro.*name=//" | perl -pe "s/ type.*$//" >interpro2name
ipr2name <- read.table("../data/interpro2name",row.names=1,col.names=c('ID','Name'),stringsAsFactors=F)

#Plot
heatdata <- list()
source("f_heatmap.2wPhylo.R")
for (i in names(clusters)) {
  cl <- clusters[[i]][['cluster']]
  #This expects a small number of clusters
  #with a bigger data set two plots must be made
  postscript(file=paste(ofold,"/pfamheat.",i,".eps",sep=''),horizontal=F,paper='special',width=9,hei=6+(nrow(restabs$pfammat)/30))
  RowSideColors=brewer.pal(max(cl),'Paired')[cl]
  col <- c("white",brewer.pal(max(data),'PuRd')) 
#  heatdata[[i]] <- heatmap.2wPhylo(x=data,Colv=tree,ColSideColors=species$Color,Rowv=cl,RowSideColors=RowSideColors,trace="none",dendrogram="column",density.info="none",margins=c(12,10),col=col,labRow=ipr2name[gsub("_.*$","",rownames(data)),])
    heatdata[[i]] <- heatmap.2wPhylo(x=data,Colv=tree,ColSideColors=species$Color,Rowv=cl,RowSideColors=RowSideColors,trace="none",dendrogram="column",density.info="none",margins=c(12,20),col=col,labRow=restabs$pfam[rownames(data),'InterproDescription'])
  dev.off()
}

#Select clustering based on pics
fclust <- 'kmeans4'
cl <- clusters[[fclust]][['cluster']]
#Join InterproPfamMAD.txt and InterproPfamOrganismStats.tx
#and clusters and principal components to output table
rows <- rownames(restabs$pfam)
outtab <- cbind(ipr2name[gsub("_.*$","",rows),],restabs$pfam,restabs$statpfam[rows,3:4],cl[rows],pca$rotation[rows,'PC1'],pca$rotation[rows,'PC2'])
colnames(outtab)[1] <- 'Name'
colnames(outtab)[grep('Median',colnames(outtab))] <- "Med of norm blast bit"
colnames(outtab)[grep('MeanAbsoluteDeviation',colnames(outtab))] <- "MAD of norm blast bit"
colnames(outtab)[grep('cl\\[rows\\]',colnames(outtab))] <- "Cluster"
colnames(outtab)[grep('PC1\"',colnames(outtab))] <- "PC1"
colnames(outtab)[grep('PC2\"',colnames(outtab))] <- "PC2"
outtab <- outtab[rev(heatdata[[fclust]][['rowInd']]),]
file <- paste(ofold,"/pfam.txt",sep='')
write.table(outtab,file=file,sep='\t',quote=F)

#Output clusters
file <- paste(ofold,"/clus.txt",sep='')
outtab <- cbind(restabs$pfam[gsub("_.*$","",rownames(restabs$clus)),'InterproDescription'],restabs$clus)
colnames(outtab)[1] <- 'InterproDescription'
corder  <- order(gsub("_.*$","",rownames(restabs$clus)),gsub("_.*$","",gsub("IPR.{6}_","" ,rownames(restabs$clus),perl=T)),as.integer(gsub("^.*_","" ,rownames(restabs$clus))))
write.table(outtab[corder,],file=file,sep='\t',quote=F)


#Special things
source("f_plotFamVsGenome.R")
postscript(file=paste(ofold,"/familyVSsize.eps",sep=''),horizontal=F,paper='special',width=6,hei=6)
reldata <- plotFamVsGenome(sp=species,pf=restabs$pfammat,highsp=c('Pestalotiopsis','Calcarisporium','Scopulariopsis'),highfam=c('IPR001138_PF00172','IPR007219_PF04082'))
dev.off()
#Table of number of clusters per species might be good too
#and that might be good for pca/clustering?

