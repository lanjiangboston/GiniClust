####################  t-test for qPCR data  ####################  
DE_t_test <- function(ExprM.RawCounts.filter,rare.cells.list.all,c_membership,out.folder,exprimentID){
  pseudo.count = 1
  data.used.log2   <- log2(ExprM.RawCounts.filter+pseudo.count)  #RNA-seq data
  #data.used.log2   <- ExprM.RawCounts.filter       #for qPCR data, because, already in log2 scale
  
  #input is a vetor in log2 space, the return is their mean in log2 space.
  Mean.in.log2space=function(x) {   
    return(log2(mean(2^(x)-1)+1))
  }
  
  for(rare.cluster in names(rare.cells.list.all) ){ 
    #rare.cluster = "Cluster_4"
    data.used3  = data.used.log2
    cells.1     = which(c_membership == rare.cluster)
    #cells.2     = which(c_membership != "cluster2")
    cells.2     = which(c_membership == "Cluster_1")
    colnames(data.used3[,cells.1])
    genes.use   = rownames(data.used3)
    mean.1      = apply(data.used3[genes.use,cells.1],1, Mean.in.log2space)
    mean.2      = apply(data.used3[genes.use,cells.2],1, Mean.in.log2space)
    total.diff  = abs(mean.1-mean.2)
    genes.diff  = names(which(total.diff> diff.cutoff))
    #genes.use   = genes.diff 
    length(genes.use) 
    p_val = unlist(lapply(genes.use,function(x)t.test(data.used3[x,cells.1],data.used3[x,cells.2])$p.value))
    log2fold_change=(mean.1-mean.2)[genes.use]
    mean_r = as.data.frame(cbind(mean.1, mean.2))
    
    differential.r = data.frame(mean.1, mean.2, p_val,log2fold_change,row.names = genes.use)
    differential.r = differential.r[with(differential.r , order(p_val, -abs(log2fold_change))), ]
    #write.table(differential.r, file=paste("results/", rare.cluster, ".diff.gene.t-test.results.xls", sep=""), sep="\t", row.names = TRUE, col.names = TRUE, quote=FALSE)
    differential.r[differential.r$p_val < 1e-22,]$p_val = 1e-22
    cluster.diffgene <- rownames(differential.r[differential.r$p_val < 1e-5 & abs(differential.r$log2fold_change) > log2(5),])
    length(cluster.diffgene) #132
    length(intersect(cluster.diffgene, GeneList.final)) #43
    area1 = length(cluster.diffgene)
    area2 = length(GeneList.final)
    n12   = length(intersect(cluster.diffgene, GeneList.final))
    cluster.overlap <- matrix(c(n12, area1-n12 , area2-n12, 283-(area1+area2-n12)),  nrow = 2)
    cluster.overlap
    fisher.test(cluster.overlap)
    
    area1 = length(cluster.diffgene)
    area2 = length(GeneList.final)
    n12   = length(intersect(cluster.diffgene, GeneList.final))
    cluster.overlap <- matrix(c(n12, area1-n12 , area2-n12, dim(ExprM.RawCounts.filter)[1]-(area1+area2-n12)),  nrow = 2)
    cluster.overlap
    ft.pvalue = fisher.test(cluster.overlap)$p.value
    
    #  VennDiagram
    if (!suppressWarnings(require("VennDiagram",quietly = TRUE))) {
      install.packages("VennDiagram", dependencies = TRUE, repos="http://cran.r-project.org")
    }
    library(VennDiagram)
    pdf(paste(out.folder, "/figures/", exprimentID,"_",rare.cluster,"_diff_gene_overlap.pdf", sep=""), height = 6, width = 6, useDingbats = FALSE)
    overlapGenes <- intersect(cluster.diffgene, GeneList.final)
    grid.newpage()
    draw.pairwise.venn(length(cluster.diffgene),
                       length(GeneList.final), 
                       length(overlapGenes),
                       category = c(paste("diffgene\n", area1), paste("HighGiniGenes\n",area2,"\n",ft.pvalue)), 
                       lty = rep("blank", 2),
                       fill = c("blue","yellow"), 
                       alpha = rep(0.5, 2), 
                       cat.pos = c(0, 0),
                       cat.dist = rep(0.025, 2))
    dev.off()
    
    #save results
    
    differential.r = differential.r[rev(order(differential.r$log2fold_change)), ]
    colnames(differential.r) = c(paste("log2.mean.", rare.cluster, sep=""), "log2.mean.Cluster_1", "t-test_pvalue", "log2fold_change" )
    write.csv(differential.r, file=paste(out.folder, "/", rare.cluster, ".diff.gene.t-test.results.csv", sep=""))
    write(overlapGenes, file=paste(out.folder, "/",rare.cluster,".overlap_genes.txt", sep=""))
  }
}
####################      End     ####################
