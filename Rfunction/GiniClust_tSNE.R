# Author: Lan Jiang
# Contact information: lan_jiang@hms.harvard.edu

#################### Visulization using tsne #################### 
GiniClust_tSNE <- function(data.type,c_membership,cell.cell.distance,perplexity,out.folder,exprimentID){
  mycols <- c("grey85",rainbow(length(unique(c_membership))-1)) 
  if(data.type == 'RNA-seq'){
    seed = 10
    set.seed(seed)  # Set a seed if you want reproducible results
    Rtsne_map <- Rtsne(cell.cell.distance, is_distance = TRUE, pca = TRUE,  max_iter = 2000, perplexity = perplexity)       
    pdf(file=paste(out.folder, "/figures/", exprimentID, "_tsne_plot.pdf", sep=""), width = 8, height =8, useDingbats = FALSE)
    plot(Rtsne_map$Y,col=mycols[as.integer(c_membership)], pch=16, xlab="Dimension_1", ylab="Dimension_2", cex=1, main="")
    legend('topleft',levels(as.factor(c_membership)),fill=mycols, bty='n', cex=1.5)
    dev.off() 
    
    # save results
    Rtnse_coord2 <-  as.data.frame(Rtsne_map$Y)
    rownames(Rtnse_coord2) = rownames(cell.cell.distance)
    colnames(Rtnse_coord2) = c("dim1", "dim2")
    write.table(Rtnse_coord2, file=paste(out.folder, "/", exprimentID,"_Rtnse_coord2.csv", sep=""), sep=",", row.names = TRUE, col.names = TRUE)
  }
  
  if(data.type == 'qPCR'){
    seed = 7
    set.seed(seed) # Set a seed if you want reproducible results
    Rtsne_map <- Rtsne(cell.cell.distance, is_distance = TRUE, pca = FALSE,  max_iter = 1000) #the result is idependent with max_iter
    pdf(file=paste(out.folder, "/figures/", exprimentID, "_tsne_plot.pdf", sep=""), width = 8, height =8, useDingbats = FALSE)
    plot(Rtsne_map$Y[,1],Rtsne_map$Y[,2],col=mycols[as.integer(c_membership)], pch=16, xlab="tSNE_1", ylab="tSNE_2", cex=1, main="")
    legend('topleft',levels(as.factor(c_membership)),fill=mycols, bty='n', cex=1.2)
    dev.off() 
    
    # save results
    Rtnse_coord2 <-  as.data.frame(Rtsne_map$Y)
    rownames(Rtnse_coord2) = rownames(cell.cell.distance)
    colnames(Rtnse_coord2) = c("dim1", "dim2")
    write.table(Rtnse_coord2, file=paste(out.folder, "/", exprimentID,"_Rtnse_coord2.csv", sep=""), sep=",", row.names = TRUE, col.names = TRUE)
  }
}
####################      End     ####################
