##############################################################################
# R source code for GiniClust v1.0 beta
# GiniClust is a computational method implemented in R to detecting rare cell types from large single cell gene expression data set. 
# It can be applied to data set come from different platform, such as multiplex qPCR data, tranditional single cell RNAseq or newly emerging UMI-based single cell RNAseq, such as inDrops and Drop-seq.
# Contact Lan Jiang (ljiang at jimmy.harvard.edu or lan.jiang.boston at gmail.com)
# Publication are comming soon. 
###############################################################################

#########################################################
### A) Installing and loading required packages
#########################################################

if (!require("ggplot2")) {
   install.packages("ggplot2", dependencies = TRUE, repos="http://cran.r-project.org")
   }
if (!require("hydroGOF")) {
   install.packages("hydroGOF", dependencies = TRUE, repos="http://cran.r-project.org")
   }
if (!require("glmnet")) {
   install.packages("glmnet", dependencies = TRUE, repos="http://cran.r-project.org")
   }
if (!require("doMC")) {
   install.packages("doMC", dependencies = TRUE, repos="http://cran.r-project.org")
   }
if (!require("StatMatch")) {
   install.packages("StatMatch", dependencies = TRUE, repos="http://cran.r-project.org")
   }
if (!require("Rtsne")) {
   install.packages("Rtsne", dependencies = TRUE, repos="http://cran.r-project.org")
   }
if (!require("fpc")) {
   install.packages("fpc", dependencies = TRUE, repos="http://cran.r-project.org")
   }
if (!require("GA")) {
   install.packages("GA", dependencies = TRUE, repos="http://cran.r-project.org")
   }
if (!require("MASS")) {
   install.packages("MASS", dependencies = TRUE, repos="http://cran.r-project.org")
   }
if (!require("session")) {
   install.packages("session", dependencies = TRUE, repos="http://cran.r-project.org")
   }
if (!require("Matrix")) {
   install.packages("Matrix", dependencies = TRUE, repos="http://cran.r-project.org")
   }
if (!require("vegan")) {
   install.packages("vegan", dependencies = TRUE, repos="http://cran.r-project.org")
   }
if (!require("data.table")) {
   install.packages("data.table", dependencies = TRUE, repos="http://cran.r-project.org")
   }
if (!require("plyr")) {
   install.packages("plyr", dependencies = TRUE, repos="http://cran.r-project.org")
   }
if (!require("reshape")) {
   install.packages("reshape", dependencies = TRUE, repos="http://cran.r-project.org")
   }
if (!require("abind")) {
   install.packages("abind", dependencies = TRUE, repos="http://cran.r-project.org")
   }
if (!require("MAST")) {
    install.packages('devtools', dependencies = TRUE, repos="http://cran.r-project.org")
    library(devtools)
    # install_github('RGLab/MAST')
    # *or* if you don't have a working latex setup
    install_github(RGLab/'MAST', build_vignettes=FALSE)
  }

#Load libraries
library(ggplot2)
library(hydroGOF)
library(glmnet)
library(doMC)
library(StatMatch)
library(Rtsne)
library(fpc)
library(GA)
library(MASS)
library(session)
library(Matrix)
library(vegan)
library(data.table)
library(plyr)
library(reshape)
library(abind)
library(MAST)


#########################################################
### B) parameters matters
#########################################################

exprimentID            = "d0t1"
CountsForNormalized    = 100000
expressed_cutoff       = 5                                   # normalized count, normal scale, not in log2 scale 
minCellNum             = 3 
minGeneNum             = 2000
max.low.cutoff         = 1                                   # filter genes based on max value, defual is 0 , in log2 scale
max.high.cutoff        = 12                                  # filter genes based on max value, defual is 12, in log2 scale 
span                   = 0.6                                 # for loess fitting             
Norm.Gini.cutoff       = 0.25                                # for choose genes which have high normlized gini score
RawCounts_cutoff       = 3
eps                    = 0.5
MinPts                 = 3
mycols                 = c("grey50", "greenyellow","red")
diff.cutoff            = 1                                   # filter gene don't have high log2_foldchange to reduce gene size in MAST analysis
lr.p_value_cutoff      = 1e-5
log2fold_change_cutoff = 1

#########################################################
### C) download data and preprocess 
#########################################################
## setwd("/groups/zhanglab/Githubupload/")
## data process before GINICLUST
#  wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1599nnn/GSM1599495/suppl/GSM1599495_ES_d0_biorep_techrep1.csv.bz2
#  bzip2 -d GSM1599495_ES_d0_biorep_techrep1.csv.bz2
#  note: files larger than 100M in this script are not uploaded to Github 

ExprM.RawCounts <- read.delim("GSM1599495_ES_d0_biorep_techrep1.csv", sep=",", head=F)    #raw data
title=c("Symbol");
for(i in 2:ncol(ExprM.RawCounts)){
  title=c(title,paste(exprimentID, ".Cell_",i-1,sep=""))
};
  colnames(ExprM.RawCounts)=title
  rownames(ExprM.RawCounts)=ExprM.RawCounts[,1]
  ExprM.RawCounts= ExprM.RawCounts[,-1]
colsum <- apply(ExprM.RawCounts,2,sum)
ExprM.normCounts <- t(t(ExprM.RawCounts)*CountsForNormalized/colsum)
write.table(ExprM.RawCounts, file=paste(exprimentID, "_rawCounts.csv", sep=""), sep=",", row.names = TRUE, col.names = TRUE, quote=FALSE)
write.table(ExprM.normCounts, file=paste(exprimentID, "_normCounts.csv", sep=""), sep=",", row.names = TRUE, col.names = TRUE, quote=FALSE)
save(ExprM.RawCounts, ExprM.normCounts, file=paste(exprimentID, "_ExprM.RData",sep=""))


#########################################################
### D) some R function
#########################################################
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

jaccard <- function(m) {
    ## common values:
    A = tcrossprod(m)
    ## indexes for non-zero common values
    im = which(A > 0, arr.ind=TRUE)
    ## counts for each row
    b = rowSums(m)

    ## only non-zero values of common
    Aim = A[im]

    ## Jacard formula: #common / (#i + #j - #common)
    J = sparseMatrix(
          i = im[,1],
          j = im[,2],
          x = Aim / (b[im[,1]] + b[im[,2]] - Aim),
          dims = dim(A)
    )

    return( J )
}

Mean.in.log2space=function(x) {   #input is a vetor in log2 space, the return is their mean in log2 space.
  return(log2(mean(2^(x)-1)+1))
}




#########################################################
### E) filter data 
#########################################################
load(paste(exprimentID, "_ExprM.RData",sep=""))
#data.used = ExprM.normCounts
ExpressedGene_per_cell=apply(ExprM.normCounts,2,function(x) length(x[x> 0 ]))
# mean(ExpressedGene_per_cell)   #[1] 5954.448
# median(ExpressedGene_per_cell) #[1] 5838
# pdf(outputfile.c)
# par(las=2,mar=c(5,3,1,1),cex=2)
# hist(ExpressedGene_per_cell,breaks=100,xlim=c(0,nrow(ExprM.normCounts)*1.1),main=paste("SingleCellData all samples\n",nrow(ExprM.normCounts)," genes vs. ",ncol(data.used)," cells\n", mean(ExpressedGene_per_cell), sep=""))
# dev.off()
ExpressedinCell_per_gene=apply(ExprM.normCounts,1,function(x) length(x[x > expressed_cutoff ]))
# mean(ExpressedinCell_per_gene)   #[1] 546.4599
# median(ExpressedinCell_per_gene) #[1] 416
# pdf(outputfile.d)
# par(las=2,mar=c(5,3,1,1),cex=2)
# hist(ExpressedinCell_per_gene,breaks=100,xlim=c(0,ncol(ExprM.normCounts)*1.1),main=paste("SingleCellData all samples\n",nrow(ExprM.normCounts)," genes vs. ",ncol(data.used)," cells\n",mean(ExpressedinCell_per_gene), sep=""))
# dev.off()
nonMir = grep("MIR|Mir", rownames(ExprM.normCounts), invert = T)  # because Mir gene is usually not accurret 
Genelist = intersect(rownames(ExprM.normCounts)[nonMir],rownames(ExprM.normCounts)[ExpressedinCell_per_gene >= minCellNum])
length(Genelist) #22830
ExprM.normCounts.filter = ExprM.normCounts[Genelist,ExpressedGene_per_cell >= minGeneNum]
ExprM.RawCounts.filter  = ExprM.RawCounts[rownames(ExprM.normCounts.filter), colnames(ExprM.normCounts.filter)]
dim(ExprM.normCounts.filter) #1] 22830  2485
dim(ExprM.RawCounts.filter)
write.table(ExprM.normCounts.filter, file=paste(exprimentID, "_gene.expression.matrix.normCounts.filtered.csv", sep=""), sep=",", row.names = TRUE, col.names = TRUE)
write.table(ExprM.RawCounts.filter,  file=paste(exprimentID, "_gene.expression.matrix.RawCounts.filtered.csv", sep=""),  sep=",", row.names = TRUE, col.names = TRUE)
save(ExprM.normCounts.filter, ExprM.RawCounts.filter, file=paste(exprimentID, "_ExprM.filter.RData", sep=""))
#22830 genes are at least expressed in 3 cell were included. Cells expressing less than 
#2,000 of these 22830 genes were excluded. The final matrix therefore contains 22830 rows (genes) quantified in 2485 samples (columns).



#########################################################
### F) generate statistics table and LOESS curve fitting 
#########################################################

# .1 load data
load(paste(exprimentID, "_ExprM.filter.RData", sep=""))

# .2 caculate (may also add more statistics in the further 10.08.2015)
Maxs        = apply(ExprM.normCounts.filter,1,max)
Gini        = apply(as.data.frame(ExprM.normCounts.filter), 1, function(x){calcul.gini(as.numeric(x)) } )    #theoretically, gini have very low chance to have a 1 value
ExprM.Stat1 = as.data.frame(cbind(Maxs,Gini))
ExprM.Stat1 = ExprM.Stat1[Maxs > 2^max.low.cutoff & Maxs < 2^max.high.cutoff , ]

# .3 fitting in max-gini space 
log2.Maxs             = log2(ExprM.Stat1$Maxs)
Gini                  = ExprM.Stat1$Gini
Gini.loess.fit        = loess(Gini~log2.Maxs, span=span, degree=1)
Normlized.Gini.Score  = residuals(Gini.loess.fit)
ExprM.Stat1           = cbind(ExprM.Stat1[,c("Maxs","Gini")], log2.Maxs, Normlized.Gini.Score)
colnames(ExprM.Stat1) = c("Maxs","Gini","log2.Maxs", "Norm.Gini")
#save(ExprM.Stat1, file= paste(exprimentID, "_ExprM.Stat1.RData", sep=""))

# .4 identify High Gini Genes 
ExprM.Stat2           = ExprM.Stat1                                        # to avoid disrupt the orignal order
ExprM.Stat2           = ExprM.Stat2[order(-ExprM.Stat2$Norm.Gini),]        # ExprM.Stat2 have beed ordered by Norm.Gini              
Genelist.HighNormGini = rownames(ExprM.Stat2[ExprM.Stat2$Norm.Gini > Norm.Gini.cutoff,])  # cutoff approach, use a cutoff to choose top genes. 

# .5 plot figures
#Histogram of gini
p1 = ggplot(ExprM.Stat1, aes(Gini, colour = "blue", fill = "blue")) 
p1 = p1 + geom_density(alpha = 0.5) 
p1 = p1 + scale_colour_manual( values = "blue" )
p1 = p1 + scale_fill_manual( values = "blue" )
p1 = p1 + labs(title="Histogram of Gini") 
p1 = p1 + theme_bw(base_size = 25, base_family = "Helvetica")                                         # clean theme
p1 = p1 + theme(panel.grid.major = element_blank(), panel.border = element_rect(colour = "black"),   
                panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
p1 = p1 + theme(legend.position="none")                                                               # remove legend
ggsave(p1, file=paste(exprimentID, "_Histogram_of_Gini.pdf", sep=""), width=8, height=8)

#Histogram of normaliezed gini
z2 = ggplot(ExprM.Stat2, aes(Norm.Gini, fill = "Norm.Gini")) + geom_density(alpha = 0.5) + labs(title="Norm.Gini") 
z2 = z2 + theme_bw(base_size = 25, base_family = "Helvetica")
z2 = z2 + theme(panel.grid.major = element_blank(), panel.border = element_rect(colour = "black"),   # clean theme
                  panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
ggsave(z2, file=paste(exprimentID, "_Histogram_of_Norm.Gini.pdf",sep=""), width=8, height=8)


#labels high normlized gini genes in log2(maximum expression value) vs gini space, smoothscatter_plot

pdf(paste(exprimentID, "_smoothscatter_plot.pdf", sep=""), width=6, height=6, useDingbats = FALSE)
smoothScatter(log2.Maxs, Gini, nrpoints = 0, xlab ="log2(Max of gene expression)", ylab = "Gini", nbin = 256, ylim=c(0,1), xlim=c(max.low.cutoff,max.high.cutoff))
points(ExprM.Stat1[Genelist.HighNormGini,]$log2.Maxs, ExprM.Stat1[Genelist.HighNormGini,]$Gini, col="red", pch=19)
#j = order(log2.Maxs)
#lines(log2.Maxs[j],Gini.loess.fit$fitted[j],col="yellow",lwd=3)
dev.off()


#barplot of high normlized gini genes with gene names, and the bar reprent normlized gini value
Stat.forbar      = ExprM.Stat2[1:length(Genelist.HighNormGini),]
Stat.forbar$Gene = factor(rownames(Stat.forbar), levels = rev(rownames(Stat.forbar))) #keep the order in the file
Stat.forbar      = as.data.frame(Stat.forbar)
b1 = ggplot(Stat.forbar, aes(x = Gene, y = Norm.Gini)) 
b1 = b1 + theme_bw(base_size = 18, base_family = "Helvetica")
b1 = b1 + theme(panel.grid.major = element_blank(), panel.border = element_rect(colour = "black"),   # clean theme
                 panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
b1 = b1 + theme(legend.position="none")   
b1 = b1 + geom_bar(stat = "identity", width=.3) + coord_flip()
ggsave(b1, file=paste(exprimentID,"_barplot_all.high.norm.gini.gene.pdf",sep=""), width=4, height=length(Genelist.HighNormGini)/3, useDingbats = FALSE)


############ .6 save results
write.table(ExprM.Stat2, file=paste(exprimentID,"_StatisticsTable_afterLOESS.csv", sep=""), sep=",", row.names = TRUE,  col.names = TRUE,  quote = FALSE) 
save(Genelist.HighNormGini, ExprM.Stat2, file=paste(exprimentID,"_StatisticsTable_afterLOESS.RData",sep=""))




#########################################################
### G) jaccard distance, dbscan for clustering
#########################################################
# .1 load data
load(paste(exprimentID,"_StatisticsTable_afterLOESS.RData",sep=""))
load(paste(exprimentID,"_ExprM.filter.RData", sep=""))

# .2 jaccard distance matrix
GeneList.final = Genelist.HighNormGini 
all.gene.as.col.ExprM.RawCounts.binary = t(ExprM.RawCounts.filter)
all.gene.as.col.ExprM.RawCounts.binary[all.gene.as.col.ExprM.RawCounts.binary <  RawCounts_cutoff] = 0
all.gene.as.col.ExprM.RawCounts.binary[all.gene.as.col.ExprM.RawCounts.binary >= RawCounts_cutoff] = 1
final.gene.as.col.ExprM.RawCounts.binary = all.gene.as.col.ExprM.RawCounts.binary[, GeneList.final]
cell.cell.jaccard.distance            = vegdist(final.gene.as.col.ExprM.RawCounts.binary, method = "jaccard")
cell.cell.jaccard.distance            = as.data.frame(as.matrix(cell.cell.jaccard.distance))
rownames(cell.cell.jaccard.distance)  = rownames(all.gene.as.col.ExprM.RawCounts.binary)
colnames(cell.cell.jaccard.distance)  = rownames(all.gene.as.col.ExprM.RawCounts.binary)

# .3 dbscan
title             = paste("eps", eps, "MinPts", MinPts, sep=".")
data.mclust       = fpc::dbscan(cell.cell.jaccard.distance, eps = eps, MinPts = MinPts,  method = "dist", showplot = FALSE)  
c_membership      = factor(paste("cluster",data.mclust$cluster ,sep=""))
levels(c_membership)[1] = "singleton"

# .4 visulization using tsne with euclidean distance
ExprM.forcluster            = t(ExprM.normCounts.filter[GeneList.final,])  #make gene as cols, cell as rows
ExprM.forcluster.log2       = log2(ExprM.forcluster+1)
euclidean_dist_matrix       = dist(ExprM.forcluster.log2, method = "euclidean")
perplexity = 30
seed = 10 
set.seed(seed) # Set a seed if you want reproducible results
Rtsne_map <- Rtsne(euclidean_dist_matrix, is_distance = TRUE, pca = FALSE,  max_iter = 1000) #the result is idependent with max_iter
pdf(file=paste(exprimentID,"_tsne_plot_euclidean_jaccard_c_membership.pdf", sep=""), width = 8, height =8, useDingbats = FALSE)
plot(Rtsne_map$Y[,1],Rtsne_map$Y[,2],col=mycols[as.integer(c_membership)], pch=16, xlab="tSNE_1", ylab="tSNE_2", cex=1, main="")
legend('topleft',levels(as.factor(c_membership)),fill=mycols, bty='n', cex=1.5)
dev.off() 

# .5 save results
Rtnse_coord2           = as.data.frame(Rtsne_map$Y)
rownames(Rtnse_coord2) = rownames(cell.cell.jaccard.distance)
colnames(Rtnse_coord2) = c("dim1", "dim2")
write.table(Rtnse_coord2, file=paste(exprimentID,"_Rtnse_coord2.csv", sep=""), sep=",", row.names = TRUE, col.names = TRUE)
save(ExprM.RawCounts.filter, cell.cell.jaccard.distance, c_membership, Rtnse_coord2, file=paste(exprimentID,"_forDE.RData",sep=""))
#note, for tSNE visualization, log2.euclidean is better than log2.scale.euclidian and log2.mahalanobis




#########################################################
### H) MAST analysis, identifing differential genes between rare cell types and other cells
#########################################################
# .1 load data
load(paste(exprimentID,"_forDE.RData", sep=""))

# .2 prepare for cluster2
data.used.df <- log2(ExprM.RawCounts.filter+1)
dim(data.used.df)       # expressioin matrix of genes(rows), cells(cols), in normal scale, not in log scale
cells.symbol.list1     = rownames(cell.cell.jaccard.distance)[which(c_membership == "cluster2")]   
cells.coord.list1      = match(cells.symbol.list1, colnames(data.used.df))                          
cells.coord.list2      = setdiff(1:dim(data.used.df)[2], cells.coord.list1)
data.1                 = data.used.df[,cells.coord.list1]
data.2                 = data.used.df[,cells.coord.list2]
mean.1                 = apply(data.1, 1, Mean.in.log2space)
mean.2                 = apply(data.2, 1, Mean.in.log2space)
mean.diff              = abs(mean.1-mean.2)
genes.list             = names(which(mean.diff> diff.cutoff))
log2fold_change        = as.data.frame(cbind(genes.list, (mean.1-mean.2)[genes.list]))
colnames(log2fold_change) = c("genename", "log2fold_change")
counts  = as.data.frame(cbind( data.1[genes.list,], data.2[genes.list,] ))
groups  = c(rep("cluster2", length(cells.coord.list1) ), rep(paste("non_", "cluster2", sep=""), length(cells.coord.list2) ) )
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
#zlm.output <- zlm.SingleCellAssay(~ Population, vbeta.1, method='glm', ebayes=TRUE)
show(zlm.output)
coefAndCI <- summary(zlm.output, logFC=TRUE)
coefAndCI <- coefAndCI[contrast != '(Intercept)',]
coefAndCI[,contrast:=abbreviate(contrast)]
zlm.lr <- lrTest(zlm.output, 'Population')
zlm.lr_pvalue <- melt(zlm.lr[,,'Pr(>Chisq)'])
zlm.lr_pvalue <- zlm.lr_pvalue[which(zlm.lr_pvalue$test.type == 'hurdle'),]
lrTest.table <-  merge(zlm.lr_pvalue, log2fold_change, by.x = "primerid", by.y = "genename")
colnames(lrTest.table) <- c("Gene", "test.type", "p_value", "log2fold_change")
cluster2_lrTest.table <- lrTest.table[order(lrTest.table$p_value),]

cluster2_two.sample <- LRT(vbeta.1, 'Population', referent='cluster2')
cluster2_MIST.r <- merge(cluster2_lrTest.table, cluster2_two.sample, by.x = "Gene", by.y = "primerid")
cluster2_MIST.r <- cluster2_MIST.r[,c("Gene","test.type.x","p_value","test.type.y","p.value", "log2fold_change")]
colnames(cluster2_MIST.r) = c("Gene","test.type.x","lr.p_value","test.type.y","LRT.p_value", "log2fold_change")
cluster2_MIST.r <- cluster2_MIST.r[order(cluster2_MIST.r$lr.p_value),]

#. 4 save results
write.table(cluster2_lrTest.table, file=paste(exprimentID,"_cluster2_lrTest.csv", sep=""), sep=",", row.names = TRUE, col.names = TRUE, quote = FALSE)
write.table(cluster2_MIST.r, file=paste(exprimentID,"_cluster2_MIST.csv", sep=""), sep=",", row.names = TRUE, col.names = TRUE, quote = FALSE)
save(cluster2_MIST.r, file=paste(exprimentID,"_cluster2_MIST.RData", sep=""))



#########################################################
### I) figures
#########################################################

# .1 vennplot plot 
load(paste(exprimentID,"_cluster2_MIST.RData", sep=""))

# .2 overlap fisher.test
rownames(cluster2_MIST.r) = cluster2_MIST.r$Gene
cluster2_MIST.r$log2fold_change = as.numeric(as.character(cluster2_MIST.r$log2fold_change)) 
cluster2.diffgene <- rownames(cluster2_MIST.r[cluster2_MIST.r$lr.p_value < lr.p_value_cutoff  & abs(cluster2_MIST.r$log2fold_change) > log2fold_change_cutoff,])
length(cluster2.diffgene)
length(intersect(cluster2.diffgene, GeneList.final)) #16
intersect(cluster2.diffgene, GeneList.final)
area1 = length(cluster2.diffgene)
area2 = length(GeneList.final)
n12   = length(intersect(cluster2.diffgene, GeneList.final))
cluster2.overlap <- matrix(c(n12, area1-n12 , area2-n12, dim(data.used.df)[1]-(area1+area2-n12)),  nrow = 2)
cluster2.overlap
fisher.test(cluster2.overlap)
# Fisher's Exact Test for Count Data
# data:  cluster2.overlap
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   43.73705 160.35596
# sample estimates:
# odds ratio 
#   85.42914 


# .3 VennDiagram
library(VennDiagram)
pdf(paste(exprimentID,"_cluster2_diff_gene_overlap.pdf", sep=""), height = 6, width = 6, useDingbats = FALSE)
overlapGenes <- intersect(cluster2.diffgene, GeneList.final)
grid.newpage()
draw.pairwise.venn(length(cluster2.diffgene),
                   length(GeneList.final), 
                   length(overlapGenes),
                   category = c("diffgene", "HighGiniGenes"), 
                   lty = rep("blank", 2),
                   fill = c("blue","yellow"), 
                   alpha = rep(0.5, 2), 
                   cat.pos = c(0, 0),
                   cat.dist = rep(0.025, 2))
dev.off()


# .4 visilization of gene expression for individual cells in orginal order .
cells.symbol.list1 = colnames(ExprM.RawCounts.filter)[which(c_membership == "cluster2")]
cells.symbol.list2 = setdiff(colnames(ExprM.RawCounts.filter),cells.symbol.list1)
for(genei in intersect(cluster2.diffgene, GeneList.final))
{
  #genei = "Gm7102";
  x1 <- as.numeric(unlist(ExprM.RawCounts.filter[genei, cells.symbol.list1]))
  x2 <- as.numeric(unlist(ExprM.RawCounts.filter[genei, cells.symbol.list2]))
  ylim = max(c(x1,x2)) * 1.2
  pdf(paste(exprimentID, "_cluster2_overlapgene_rawCounts_bar_plot.", genei, ".pdf", sep=""), height = 6, width = 12, useDingbats = FALSE)
  par(mfrow=c(1,2))
  barplot(x1, main = paste(genei, "cluster2", sep=" "), ylim=c(0,ylim), col=mycols[3], xlab="", ylab="UMI/cell", width = 0.1, axisnames = FALSE, space = 1.9);
  barplot(x2, main = "others", ylim=c(0,ylim), width = 1, col=mycols[2], xlab="", axisnames = FALSE);
  dev.off()

}


table(c_membership)
# c_membership
#    singleton cluster1 cluster2 
#        1     2480        4 

library("session")
save.session(paste(exprimentID, "_", format(Sys.time(), "%b%d"),".session",sep=""))
#restore.session("results/Dec")





