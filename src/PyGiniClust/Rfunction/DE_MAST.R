# Author: Lan Jiang
# Contact information: lan_jiang@hms.harvard.edu

####################  MAST analysis for RNA seq ####################
DE_MAST <- function(ExprM.RawCounts.filter,rare.cells.list.all,out.folder,exprimentID){
  #identifing differential genes between rare cell types and other cells
  
  #for RNA-seq only consider Rare cell type up-regulated genes, because we use one-direction gini.
  #' input is a data matrix, which is already in log2 scale, 
  #' and group vector, like c(1,1,1,0,0), 
  #' and pseudo.count, which is used for log2 transformation
  
  #' output is three column, log2(mean.1), log2(mean.2), log2_foldchange.
  ####### function data log2 mean, and log2 fold changle ######
  #install.packages("varhandle")
  if (!suppressWarnings(require("MAST",quietly = TRUE))) {
    if(!suppressWarnings(require("RCurl",quietly = TRUE))) {
      install.packages('RCurl', dependencies = TRUE, repos="http://cran.r-project.org")
    }
    library(RCurl)
    library(httr)
    set_config( config( ssl_verifypeer = 0L ) )
    if(!suppressWarnings(require("devtools",quietly = TRUE))) {
      install.packages('devtools', dependencies = TRUE, repos="http://cran.r-project.org")
    }
    library(devtools)
    source("http://bioconductor.org/biocLite.R")
    biocLite("Biobase",suppressUpdates=TRUE,suppressAutoUpdate=TRUE)
    biocLite("BiocGenerics",suppressUpdates=TRUE,suppressAutoUpdate=TRUE)
    
    # install_github('RGLab/MAST')
    # *or* if you don't have a working latex setup
    install_github('RGLab/MAST', build_vignettes=FALSE)
  }
  library(MAST)
  
  if (!suppressWarnings(require("varhandle",quietly = TRUE))) {
    install.packages("varhandle", dependencies = TRUE, repos="http://cran.r-project.org")
  }
  library(varhandle)
  Mean.in.log2space=function(x,pseudo.count) {
    return(log2(mean(2^(x)-pseudo.count)+pseudo.count))
  }
  
  stat.log2=function(data.m, group.v, pseudo.count){
    #data.m=data.used.log2
    log2.mean.r <- aggregate(t(data.m), list(as.character(group.v)), function(x) Mean.in.log2space(x,pseudo.count))
    log2.mean.r <- t(log2.mean.r)
    colnames(log2.mean.r) <- paste("mean.group",log2.mean.r[1,], sep="")
    log2.mean.r = log2.mean.r[-1,]
    log2.mean.r = as.data.frame(log2.mean.r)
    log2.mean.r = unfactor(log2.mean.r)  #from varhandle
    log2.mean.r[,1] = as.numeric(log2.mean.r[,1])
    log2.mean.r[,2] = as.numeric(log2.mean.r[,2])
    log2_foldchange = log2.mean.r$mean.group1-log2.mean.r$mean.group0
    results = data.frame(cbind(log2.mean.r$mean.group0,log2.mean.r$mean.group1,log2_foldchange))
    colnames(results) = c("log2.mean.group0","log2.mean.group1","log2_fc")
    rownames(results) = rownames(log2.mean.r)
    return(results)
  }
  
  ####### function END ######
  
  
  ####### function m.auc  ######
  #install.packages("ROCR")
  if (!suppressWarnings(require("ROCR",quietly = TRUE))) {
    install.packages("ROCR", dependencies = TRUE, repos="http://cran.r-project.org")
  }
  library("ROCR")
  v.auc = function(data.v,group.v) {
    prediction.use=prediction(data.v, group.v, 0:1)
    perf.use=performance(prediction.use,"auc")
    auc.use=round(perf.use@y.values[[1]],3)
    return(auc.use)
  }
  m.auc=function(data.m,group.v) {
    AUC=apply(data.m, 1, function(x) v.auc(x,group.v))
    AUC[is.na(AUC)]=0.5
    return(AUC)
    
  }  
  ####### function m.auc END ######
  
  
  #input
  pseudo.count = 0.1
  data.used.log2   <- log2(ExprM.RawCounts.filter+pseudo.count)  #RNA-seq data
  
  ##loop for all rare cell type
  for(rare.cluster in names(rare.cells.list.all) ){
    
    cells.symbol.list2     = rare.cells.list.all[[rare.cluster]]
    cells.coord.list2      = match(cells.symbol.list2, colnames(data.used.log2))                          
    cells.symbol.list1     = clustering_membership_r$cell.ID[which(clustering_membership_r$cluster.ID == "Cluster_1")];  #always compare with the biggest cluster, after rename
    cells.coord.list1      = match(cells.symbol.list1, colnames(data.used.log2))   
    data.used.log2.ordered  = cbind(data.used.log2[,cells.coord.list1], data.used.log2[,cells.coord.list2])
    group.v <- c(rep(0,length(cells.coord.list1)), rep(1, length(cells.coord.list2)))
    #ouput
    log2.stat.result <- stat.log2(data.used.log2.ordered, group.v, pseudo.count)
    Auc <- m.auc(data.used.log2.ordered, group.v)
    bigtable <- data.frame(cbind(log2.stat.result, Auc))
    
    diff.cutoff = 1
    DE <- bigtable[bigtable$log2_fc >diff.cutoff,] 
    dim(DE)
    data.1                 = data.used.log2[,cells.coord.list1]
    data.2                 = data.used.log2[,cells.coord.list2]
    genes.list = rownames(DE)
    log2fold_change        = cbind(genes.list, DE$log2_fc)
    colnames(log2fold_change) = c("gene.name", "log2fold_change")
    counts  = as.data.frame(cbind( data.1[genes.list,], data.2[genes.list,] ))
    #groups  = c(rep(rare.cluster, length(cells.coord.list1) ), rep(paste("non_", rare.cluster, sep=""), length(cells.coord.list2) ) )
    groups  = c(rep("Cluster_1", length(cells.coord.list1) ), rep(rare.cluster, length(cells.coord.list2) ) )
    groups  = as.character(groups)
    data_for_MIST <- as.data.frame(cbind(rep(rownames(counts), dim(counts)[2]), melt(counts),rep(groups, each = dim(counts)[1]), rep(1, dim(counts)[1] * dim(counts)[2]) ))
    colnames(data_for_MIST) = c("Gene", "Subject.ID", "Et", "Population", "Number.of.Cells")
    vbeta = data_for_MIST
    vbeta.fa <- FluidigmAssay(vbeta, idvars=c("Subject.ID"),
                              primerid='Gene', measurement='Et', ncells='Number.of.Cells',
                              geneid="Gene",  cellvars=c('Number.of.Cells', 'Population'),
                              phenovars=c('Population'), id='vbeta all')
    vbeta.1 <- subset(vbeta.fa, ncells==1)
    # .3 MAST 
    layername(vbeta.1)
    head(cData(vbeta.1))
    zlm.output <- zlm.SingleCellAssay(~ Population, vbeta.1, method='bayesglm', ebayes=TRUE)
    show(zlm.output)
    coefAndCI <- summary(zlm.output, logFC=TRUE)
    coefAndCI <- coefAndCI[contrast != '(Intercept)',]
    coefAndCI[,contrast:=abbreviate(contrast)]
    zlm.lr <- lrTest(zlm.output, 'Population')
    zlm.lr_pvalue <- melt(zlm.lr[,,'Pr(>Chisq)'])
    zlm.lr_pvalue <- zlm.lr_pvalue[which(zlm.lr_pvalue$test.type == 'hurdle'),]
    
    
    
    lrTest.table <-  merge(zlm.lr_pvalue, DE, by.x = "primerid", by.y = "row.names")
    colnames(lrTest.table) <- c("Gene", "test.type", "p_value", paste("log2.mean.", "Cluster_1", sep=""), paste("log2.mean.",rare.cluster,sep=""), "log2fold_change", "Auc")
    cluster_lrTest.table <- lrTest.table[rev(order(lrTest.table$Auc)),]
    
    #. 4 save results
    write.csv(cluster_lrTest.table, file=paste(out.folder, "/",rare.cluster,".diff.gene.lrTest.results.csv", sep=""))
  }
  return(cluster_lrTest.table)
}
####################      End     ####################
