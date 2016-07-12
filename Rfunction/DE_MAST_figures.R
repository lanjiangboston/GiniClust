####################  MAST figures ####################      
DE_MAST_figures <- function(rare.cells.list.all,cluster_lrTest.table,clustering_membership_r,out.folder,exprimentID){
  for(rare.cluster in names(rare.cells.list.all)){
    # .2 overlap fisher.test
    rownames(cluster_lrTest.table) = cluster_lrTest.table$Gene
    cluster_lrTest.table$log2fold_change = as.numeric(as.character(cluster_lrTest.table$log2fold_change)) 
    cluster.diffgene <- rownames(cluster_lrTest.table[cluster_lrTest.table$p_value < lr.p_value_cutoff  & abs(cluster_lrTest.table$log2fold_change) > diff.cutoff,])
    length(cluster.diffgene)
    length(intersect(cluster.diffgene, GeneList.final)) 
    intersect(cluster.diffgene, GeneList.final)
    
    area1 = length(cluster.diffgene)
    area2 = length(GeneList.final)
    n12   = length(intersect(cluster.diffgene, GeneList.final))
    cluster.overlap <- matrix(c(n12, area1-n12 , area2-n12, dim(ExprM.RawCounts.filter)[1]-(area1+area2-n12)),  nrow = 2)
    cluster.overlap
    ft.pvalue = fisher.test(cluster.overlap)$p.value
    
    # .3 VennDiagram
    if (!suppressWarnings(require("VennDiagram",quietly = TRUE))) {
      install.packages("VennDiagram", dependencies = TRUE, repos="http://cran.r-project.org")
    }
    library(VennDiagram)
    pdf(paste(out.folder,"/figures/", exprimentID,"_",rare.cluster,"_diff_gene_overlap.pdf", sep=""), height = 6, width = 6, useDingbats = FALSE)
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
    
    
    # .4 barplot visilization of gene expression for individual cells in orginal order .
    cells.symbol.list1 = rare.cells.list.all[[rare.cluster]]
    cells.symbol.list2     = clustering_membership_r$cell.ID[which(clustering_membership_r$cluster.ID == "Cluster_1")];  # compare with the biggest cluster, after rename 
    genelist.forbarplot = intersect( rownames(cluster_lrTest.table[cluster_lrTest.table$Auc>0.98,]), GeneList.final )
    
    if(length(overlapGenes) > 0)
    {
      for(genei in overlapGenes)
      {
        x1 <- as.numeric(unlist(ExprM.RawCounts.filter[genei, cells.symbol.list1]))
        x2 <- as.numeric(unlist(ExprM.RawCounts.filter[genei, cells.symbol.list2]))
        ylim = max(c(x1,x2)) * 1.2
        pdf(paste(out.folder, "/figures/", exprimentID, "_",rare.cluster, "_overlapgene_rawCounts_bar_plot.", genei, ".pdf", sep=""), height = 6, width = 12, useDingbats = FALSE)
        par(mfrow=c(1,2))
        barplot(x1, main = paste(genei, rare.cluster, sep=" "), ylim=c(0,ylim), col=mycols[3], xlab="", ylab="Counts/cell", width = 0.1, axisnames = FALSE, space = 1.9);
        barplot(x2, main = "Cluster_1", ylim=c(0,ylim), width = 1, col=mycols[2], xlab="", axisnames = FALSE);
        dev.off()
      }
    }
  }
}
####################      End     ####################
