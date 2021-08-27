args=(commandArgs(TRUE))

#Read counts table from command line
ct <- args[1]
ct <- read.table(ct, sep=",", header=T, row.names=1, stringsAsFactors = F, check.names = F)
mu <- args[2]


#Create Dendrogram
b_frame<-as.data.frame(ct)
new_frame<-apply(b_frame, 2, as.numeric)
colnames(new_frame)<-substr(colnames(new_frame),1,8) 
dat<-as.matrix(new_frame) 
cd<-dist(t(dat)) 
cc<-hclust(cd) 
hcd = as.dendrogram(cc) 

pdf("dendrogram.pdf")
par(cex=0.8) 
plot(hcd) 
dev.off()  
cat("Dendrogram Completed\n")


#If mapped_unmapped table is present - create bar plot
if(!is.null(mu)){
  mu <- read.table(mu, sep=",", header=T, row.names=1, stringsAsFactors = F, check.names = F)
  
  #Generate bar plot of total reads/sample and the mapped and unmapped portion.
  ll2<-as.matrix(mu)
  
  mettol<-seq(1,dim(ll2)[2],50)
  pdf("barplot_mapped_unmapped.pdf")
  par(cex=1.5) 
  par(las=2)
  
  for (i in 1:length(mettol)){
    if(dim(ll2)[2]>=(mettol[i]+49)){
      barplot(as.matrix(ll2[,mettol[i]:(mettol[i]+49)]),  col=c("yellow", "orange"), ylab="read counts",xlab="sample name", cex.axis=0.5, cex.names=0.5) 
      legend("topleft", c("unmapped","mapped"), fill=c("yellow", "orange"), cex=0.7)
    }else{
      barplot(as.matrix(ll2[,mettol[i]:dim(ll2)[2]]),  col=c("yellow", "orange"), ylab="read counts",xlab="sample name", cex.axis=0.5, cex.names=0.5) 
      legend("topleft", c("unmapped","mapped"), fill=c("yellow", "orange"), cex=0.7) 
    }
  }
  dev.off()
  cat("Bar plot Completed\n")
  
}


#Create PCA
require(DESeq2)
require(ggplot2)


#Row Variance - taken from genefilter package
rowVars <- function (x, ...){
  sqr = function(x) x * x
  n = rowSums(!is.na(x))
  n[n <= 1] = NA
  return(rowSums(sqr(x - rowMeans(x, ...)), ...)/(n - 1))
}

#Round up to nearest 10, 100, etc
RoundUp <- function(from,to) ceiling(from/to)*to


#Variance Stabilization Function - from DESeq package 
varStab <- function (object, blind = TRUE, nsub, fitType = "parametric") {
  if (nrow(object) < nsub) {
    stop("less than 'nsub' rows,\n  it is recommended to use varianceStabilizingTransformation directly")
  }
  if (is.null(colnames(object))) {
    colnames(object) <- seq_len(ncol(object))
  }
  if (is.matrix(object)) {
    matrixIn <- TRUE
    object <- DESeqDataSetFromMatrix(object, DataFrame(row.names = colnames(object)), 
                                     ~1)
  }
  else {
    if (blind) {
      design(object) <- ~1
    }
    matrixIn <- FALSE
  }
  if (is.null(sizeFactors(object)) & is.null(normalizationFactors(object))) {
    object <- estimateSizeFactors(object)
  }
  baseMean <- rowMeans(counts(object, normalized = TRUE))
  if (sum(baseMean > 5) < nsub) {
    stop("less than 'nsub' rows with mean normalized count > 5, \n  it is recommended to use varianceStabilizingTransformation directly")
  }
  object.sub <- object[baseMean > 5, ]
  baseMean <- baseMean[baseMean > 5]
  o <- order(baseMean)
  idx <- o[round(seq(from = 1, to = length(o), length = nsub))]
  object.sub <- object.sub[idx, ]
  object.sub <- estimateDispersionsGeneEst(object.sub, quiet = TRUE)
  object.sub <- estimateDispersionsFit(object.sub, fitType = fitType, 
                                       quiet = TRUE)
  suppressMessages({
    dispersionFunction(object) <- dispersionFunction(object.sub)
  })
  vsd <- varianceStabilizingTransformation(object, blind = FALSE)
  if (matrixIn) {
    return(assay(vsd))
  }
  else {
    return(vsd)
  }
}


#PCA pre-process function
pca_fun <- function(x, condition){
  #create column data
  coldata <- data.frame(row.names=colnames(x), condition)
  #make DESeqDataSetFromMatrix
  dds<-DESeqDataSetFromMatrix(countData=b_frame, colData=coldata, design=~condition)
  return(dds)
}


#process data for PCA plot
condition <-factor(c(colnames(ct)))
f_dds <- pca_fun(ct, condition)
ntop <- 500
nsub <- RoundUp(ncol(ct)*0.90, 1)
vsd <- varStab(f_dds, nsub=nsub)
rv <- rowVars(assay(vsd))
select = order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
pca = prcomp(t(assay(vsd)[select, ]))
scores <- data.frame(condition, pca$x)


#Plot PCA and output to PDF
pca_input <- as.vector(colnames(ct))
pca_in <- strsplit(pca_input, ":")
pdf("PCA.pdf")
#Plot PC1 vs. PC2
#To plot different principal components, choose different columns from dataframe 'scores' 
pca_p <- qplot(x=scores[,2], y=scores[,3], data=scores, 
               colour=factor(condition), xlab="PC1", ylab="PC2", size=I(3.0)) + labs(colour='Samples')
print(pca_p)
dev.off()
cat("PCA Completed\n")



