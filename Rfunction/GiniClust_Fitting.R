# Author: Lan Jiang
# Contact information: lan_jiang@hms.harvard.edu

#################### LOESS curve fitting #################### 
GiniClust_Fitting <- function(data.type,ExprM.RawCounts.filter,out.folder,exprimentID){
  
  # Calculate gini index
  calcul.gini = function(x, unbiased = TRUE, na.rm = FALSE){
    if (!is.numeric(x)){
      warning("'x' is not numeric; returning NA")
      return(NA)
    }
    if (!na.rm && any(na.ind = is.na(x)))
      stop("'x' contain NAs")
    if (na.rm)
      x = x[!na.ind]
    n = length(x)
    mu = mean(x)
    N = if (unbiased) n * (n - 1) else n * n
    ox = x[order(x)]
    dsum = drop(crossprod(2 * 1:n - n - 1,  ox))
    dsum / (mu * N)
  }
  
  if(data.type == 'RNA-seq'){
    gini = apply(as.data.frame(ExprM.RawCounts.filter), 1, function(x){calcul.gini(as.numeric(x)) } )    #theoretically, gini have very low chance to have a 1 value
    GiniIndex = as.data.frame(cbind(1:dim(ExprM.RawCounts.filter)[1], gini))
  } 
  if(data.type == 'qPCR'){
    GiniIndex1 <- as.data.frame(apply(ExprM.RawCounts.filter, 1, function(x){calcul.gini(as.numeric(x)) } ) )
    GiniIndex2 <- as.data.frame(apply(ExprM.RawCounts.filter+0.00001, 1, function(x){calcul.gini(as.numeric(1/x)) } ) ) #bi directional
    GiniIndex  <- cbind(GiniIndex1, GiniIndex2)
    colnames(GiniIndex)=c("gini1","gini2")
    GiniIndex$gini2_sign = 0 - GiniIndex$gini2;
    GiniIndex$gini = apply(GiniIndex, 1, max)
    GiniIndex <- na.omit(GiniIndex)
    GiniIndex$gini_sign = GiniIndex$gini
    for(genei in 1:dim(GiniIndex)[1])
    {
      GiniIndex[genei, 5] = ifelse(  GiniIndex[genei, 1] > GiniIndex[genei,2], "up-regulation", "down-regulation") 
    }
    #dim(GiniIndex) 
    write.table(GiniIndex, file=paste(out.folder, "/", exprimentID,"_bi-directional.GiniIndexTable.csv", sep=""), sep=",", row.names = TRUE,  col.names = TRUE,  quote = FALSE) 
  }
  
  Maxs          = apply(ExprM.RawCounts.filter,1,max)
  Means         = apply(ExprM.RawCounts.filter,1,mean)
  log2.Maxs     = log2(Maxs+0.1)
  ExprM.Stat1   = as.data.frame(cbind(Maxs,GiniIndex$gini,log2.Maxs))
  colnames(ExprM.Stat1) = c("Maxs","Gini","log2.Maxs")
  ExprM.Stat1 = ExprM.Stat1[ExprM.Stat1$log2.Maxs>log2.expr.cutoffl & ExprM.Stat1$log2.Maxs<=log2.expr.cutoffh ,] 
  log2.Maxs = ExprM.Stat1$log2.Maxs
  Gini      = ExprM.Stat1$Gini
  Maxs      = ExprM.Stat1$Maxs
  
  # fitting in max-gini space 
  Gini.loess.fit        = loess(Gini~log2.Maxs, span=span, degree=1)
  Normlized.Gini.Score  = Gini.loess.fit$residuals   #residuals = Gini - Gini.fitted
  Gini.fitted           = Gini.loess.fit$fitted    
  ExprM.Stat1           = as.data.frame(cbind(ExprM.Stat1[,c("Maxs","Gini", "log2.Maxs")], Normlized.Gini.Score, Gini.fitted))
  colnames(ExprM.Stat1) = c("Maxs","Gini","log2.Maxs", "Norm.Gini", "Gini.fitted")
  ### remove 25% of first round outlier genes, do second round loess
  Gini.loess.fit.residual = residuals(Gini.loess.fit)                               
  thresh.outlier = quantile(Gini.loess.fit.residual[Gini.loess.fit.residual>0], outlier_remove) 
  id.genes.loess.fit = which(Gini.loess.fit.residual < thresh.outlier)               
  id.outliers.loess.fit = which(Gini.loess.fit.residual >= thresh.outlier)          
  log2.Maxs.genes = log2.Maxs[id.genes.loess.fit]                                   
  log2.Maxs.outliers = log2.Maxs[id.outliers.loess.fit]                            
  Gini.loess.fit.2 = loess(Gini[id.genes.loess.fit]~log2.Maxs[id.genes.loess.fit], span=span, degree = 1)
  Gini.loess.fit.2.predict = predict(Gini.loess.fit.2)  
  
  #plot second round fit
  Gini.loess.fit.2.x.y = cbind(log2.Maxs.genes,Gini.loess.fit.2.predict)
  Gini.loess.fit.2.x.y.uniq = as.data.frame(unique(Gini.loess.fit.2.x.y))
  Gini.loess.fit.2.x.y.uniq = Gini.loess.fit.2.x.y.uniq[order(Gini.loess.fit.2.x.y.uniq[,1]),]
  log2.Maxs.genes.sorted = log2.Maxs.genes[order(log2.Maxs.genes)]                   
  Gini.loess.fit.2.predict.sorted = Gini.loess.fit.2.predict[order(log2.Maxs.genes)] 
  #using Gini.loess.fit.2 as model, predict gini value for those outlier which are not used for build model.
  #for each max in outliers set, find the id of max value which is most close in fitted data set
  loc.outliers = apply(matrix(log2.Maxs.outliers),1,function(x){
    if(x<max(log2.Maxs.genes.sorted)){
      return(which(log2.Maxs.genes.sorted>=x)[1])
    }else{
      return(which.max(log2.Maxs.genes.sorted))
    }})                
  #check the results
  outlier_max_in_fit <- cbind(log2.Maxs.outliers, loc.outliers, log2.Maxs.genes.sorted[loc.outliers])
  
  #based on Gini.loess.fit.2, predict outliers which was not used for fitting
  Gini.outliers.predict = apply(cbind(seq(length(log2.Maxs.outliers)),log2.Maxs.outliers),1,function(x){
    id = x[1]
    value = x[2]
    if(value == log2.Maxs.genes.sorted[loc.outliers[id]]){
      return(as.numeric(Gini.loess.fit.2.x.y.uniq[which(Gini.loess.fit.2.x.y.uniq$log2.Maxs.genes>=value)[1],2]))
    }else{
      if(loc.outliers[id]>1){
        return(Gini.loess.fit.2.predict.sorted[loc.outliers[id]-1]+(Gini.loess.fit.2.predict.sorted[loc.outliers[id]]-Gini.loess.fit.2.predict.sorted[loc.outliers[id]-1])*(value-log2.Maxs.genes.sorted[loc.outliers[id]-1])/(log2.Maxs.genes.sorted[loc.outliers[id]]-log2.Maxs.genes.sorted[loc.outliers[id]-1]))
      }else{
        return(Gini.loess.fit.2.predict.sorted[2]-(Gini.loess.fit.2.predict.sorted[2]-Gini.loess.fit.2.predict.sorted[1])*(log2.Maxs.genes.sorted[2]-value)/(log2.Maxs.genes.sorted[2]-log2.Maxs.genes.sorted[1]))
      }
    }
  })
  
  #plot outliers predict results
  outliers.precit.x.y.uniq = as.data.frame(unique(cbind(log2.Maxs.outliers, Gini.outliers.predict)))
  #plot(outliers.precit.x.y.uniq)
  #plot whole fit2 
  colnames(outliers.precit.x.y.uniq) = colnames(Gini.loess.fit.2.x.y.uniq)
  Gini.loess.fit.2.full.x.y.uniq = rbind(Gini.loess.fit.2.x.y.uniq, outliers.precit.x.y.uniq)
  #plot(Gini.loess.fit.2.full.x.y.uniq)
  
  #calcualte Normlized.Gini.Score2
  Normlized.Gini.Score2                        = rep(0,length(Gini.loess.fit.residual))               
  Normlized.Gini.Score2[id.genes.loess.fit]    = residuals(Gini.loess.fit.2)                         
  Normlized.Gini.Score2[id.outliers.loess.fit] = Gini[id.outliers.loess.fit] - Gini.outliers.predict 
  
  Gini.fitted2           = Gini - Normlized.Gini.Score2         
  ExprM.Stat1            = as.data.frame(cbind(ExprM.Stat1[,c("Maxs","Gini", "log2.Maxs", "Gini.fitted", "Norm.Gini" )], Gini.fitted2, Normlized.Gini.Score2))
  colnames(ExprM.Stat1)  = c("Maxs","Gini","log2.Maxs", "Gini.fitted","Norm.Gini",  "Gini.fitted2", "Norm.Gini2")
  Gini.pvalue            = pnorm(-abs(scale(ExprM.Stat1$Norm.Gini2, center=TRUE,scale=TRUE)))
  ExprM.Stat2            = cbind(ExprM.Stat1, Gini.pvalue)  #first time use ExprM.Stat2
  
  #for each measurement, first ranked by themself.
  # dentify High Gini Genes with Norm.Gini
  ExprM.Stat2         = ExprM.Stat2[rev(order(ExprM.Stat2$Norm.Gini2)),]
  Genelist.HighNormGini = rownames(ExprM.Stat2[ExprM.Stat2$Norm.Gini2 > Norm.Gini.cutoff,])  # cutoff approach, use a cutoff to choose top genes. 
  length(Genelist.HighNormGini)
  
  # identify High Gini Genes with pvalue
  ExprM.Stat2         = ExprM.Stat2[order(ExprM.Stat2$Gini.pvalue),]
  Genelist.top_pvalue = rownames(ExprM.Stat2[ExprM.Stat2$Gini.pvalue < Gini.pvalue_cutoff & ExprM.Stat2$Norm.Gini2 > 0,])
  length(Genelist.top_pvalue)
  
  
  # plot figures
  xall = ExprM.Stat2$log2.Maxs
  yall = ExprM.Stat2$Gini
  yfit2 = Gini.fitted2 
  #Histogram of pvalue
  main=paste("Histogram of -log10(Gini.pvalue)", "\ncutoff=",Gini.pvalue_cutoff, "\n Gene num=", length(Genelist.top_pvalue), sep="")
  pdf(paste(out.folder, "/figures/", exprimentID, "_histogram of Normalized.Gini.Socre.pdf", sep=""), width=6, height=6)
  hist(0-log10(ExprM.Stat2$Gini.pvalue), breaks=100, main=main)
  abline(v=0-log10(Gini.pvalue_cutoff), col="red")
  dev.off()
  
  # save results
  
  if(dim(GiniIndex)[2] > 2){
    Gini_related_table <- merge(GiniIndex, ExprM.Stat2, by.x = "row.names", by.y = "row.names")
    Gini_related_table_qPCR <- Gini_related_table[,c("Row.names", "log2.Maxs", "Gini", "gini_sign", "Norm.Gini2")]
    #Sign of Gini index. 1 for up-regulation and -1 for down-regulation.
    Gini_related_table_qPCR = Gini_related_table_qPCR[rev(order(Gini_related_table_qPCR$Norm.Gini2)),]
    colnames(Gini_related_table_qPCR) = c("Gene","log2( Max Expression Level )", "Raw Gini Index", "Sign of Gini index", "Normalized Gini Index")
    write.csv(Gini_related_table_qPCR, file=paste(out.folder,"/Gini_related_table_qPCR.csv",sep=""))
    
    pdf(paste(out.folder, "/figures/", exprimentID, "_smoothScatter_normGini_genes.pdf", sep=""), width=6, height=6, useDingbats = FALSE)
    smoothScatter(xall, yall, nrpoints = 0, xlab ="log2( Max Expression Level )", ylab = "Gini", main=paste("Gene num=",length(Genelist.HighNormGini),"\nHighNormGini=",Norm.Gini.cutoff,sep=""), nbin = 256, ylim=c(0,1), xlim=c(0,max(xall)))
    points(ExprM.Stat2[Genelist.HighNormGini,]$log2.Maxs, ExprM.Stat2[Genelist.HighNormGini,]$Gini, col="red", pch=19)
    dev.off()
    
  }else{
    Gini_related_table <- ExprM.Stat2
    Gini_related_table_RNAseq <- Gini_related_table[,c( "log2.Maxs","Gini", "Norm.Gini2", "Gini.pvalue")]
    Gini_related_table_RNAseq <- cbind(rownames(Gini_related_table_RNAseq ), Gini_related_table_RNAseq )
    colnames(Gini_related_table_RNAseq) <- c("Gene", "log2( Max Expression Level )", "Raw Gini Index", "Normalized Gini Index", "p-value")
    write.csv(Gini_related_table_RNAseq, file=paste(out.folder,"/Gini_related_table_RNA-seq.csv",sep=""))
    
    pdf(paste(out.folder, "/figures/", exprimentID, "_smoothScatter_pvalue_gene.pdf", sep=""), width=6, height=6, useDingbats = FALSE)
    smoothScatter(xall, yall, nrpoints = 0, xlab ="log2( Max Expression Level )", ylab = "Gini", main=paste("Gene num=",length(Genelist.top_pvalue),"\npvalue=",Gini.pvalue_cutoff,sep=""), nbin = 256, ylim=c(0,1), xlim=c(0,max(xall)))
    points(ExprM.Stat2[Genelist.top_pvalue,]$log2.Maxs, ExprM.Stat2[Genelist.top_pvalue,]$Gini, col="red", pch=19)
    dev.off()
  }
  if(data.type == 'RNA-seq'){
    return(Genelist.top_pvalue)
  }
  if(data.type == 'qPCR'){
    return(Genelist.HighNormGini)
  }
}
####################      End     ####################
