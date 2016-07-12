#################### dbscan for clustering ####################
GiniClust_Clustering <- function(data.type,ExprM.RawCounts.filter,GeneList.final,eps,MinPts,out.folder,exprimentID){
  
  # Calculate jaccard distances
  jaccard <- function(m){
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
  
  if(data.type == 'RNA-seq'){
    # binarization
    m = ExprM.RawCounts.filter[GeneList.final,]
    m2 = m
    # if t is the expression vector, the question asked here is which count value x is the smallest one, when sum(t[t>x]) > Gamma % * sum(t). 
    bc.list.low = c()
    bc.list.med = c()
    bc.list.high = c()
    for(rn in 1:dim(m)[1])
    {
      t <- as.numeric(m[rn,])
      t.table <- data.frame(table(t))
      c=as.numeric(as.character(t.table$t))
      f=as.numeric(as.character(t.table$Freq))
      tcf  <- data.frame(cbind(c,f))
      csum <- apply(tcf,1,function(x){sum(t[t>=x[1]])/sum(t[t>0])})
      tcfs <- data.frame(cbind(c,f,csum))
      tcfs <- tcfs[rev(order(tcfs$c)),]
      n = max(3, which(tcfs$csum>Gamma)[1])
      bc.list.high[rn]  = tcfs$c[n]
      bc.list.low[rn] = tcfs$c[n+1]   #range will be > bc.list.low[rn] to bc.list.high[rn] 
      bc.list.med[rn] = mean(c(bc.list.high[rn], bc.list.low[rn]))
    }
    top_n_gene = max(length(bc.list.low)*0.1, 10)
    RawCounts_cutoff = floor(mean(bc.list.med[1:top_n_gene]))
    #binarization
    all.gene.as.col.ExprM.RawCounts.binary = t(ExprM.RawCounts.filter)
    all.gene.as.col.ExprM.RawCounts.binary[all.gene.as.col.ExprM.RawCounts.binary <  RawCounts_cutoff] = 0
    all.gene.as.col.ExprM.RawCounts.binary[all.gene.as.col.ExprM.RawCounts.binary >= RawCounts_cutoff] = 1
    final.gene.as.col.ExprM.RawCounts.binary = all.gene.as.col.ExprM.RawCounts.binary[, intersect(colnames(all.gene.as.col.ExprM.RawCounts.binary), GeneList.final)]
    dim(final.gene.as.col.ExprM.RawCounts.binary) 
    
    # jaccard distance matrix
    #locate cells whose gene expression values are all zeros
    index.cell.zero = which(apply(final.gene.as.col.ExprM.RawCounts.binary,1,function(x) length(which(x>0)))==0)
    cell.cell.jaccard.distance = 1 - jaccard(final.gene.as.col.ExprM.RawCounts.binary)
    cell.cell.jaccard.distance = as.data.frame(as.matrix(cell.cell.jaccard.distance))
    #covert distance between two 'zero cells' to zero
    cell.cell.jaccard.distance[index.cell.zero,index.cell.zero] = 0
    cell.cell.jaccard.distance            = as.data.frame(as.matrix(cell.cell.jaccard.distance))
    rownames(cell.cell.jaccard.distance)  = rownames(final.gene.as.col.ExprM.RawCounts.binary)
    colnames(cell.cell.jaccard.distance)  = rownames(final.gene.as.col.ExprM.RawCounts.binary)
    dim(cell.cell.jaccard.distance)
    ###-----------------------------------------------------------------------------------###
    
    
    # dbscan
    title             = paste("eps", eps, "MinPts", MinPts, sep=".")
    data.mclust       = fpc::dbscan(cell.cell.jaccard.distance, eps = eps, MinPts = MinPts,  method = "dist", showplot = FALSE)  
    
    #rename cluster names based on the size of each cluster.
    o_membership      = factor(paste("db_",data.mclust$cluster ,sep=""))
    if(levels(o_membership)[1]=="db_0"){levels(o_membership)[1] = "Singleton"}
    c_membership = o_membership
    cluster_stat = as.data.frame(table(o_membership))
    cluster_stat = cluster_stat[cluster_stat$o_membership != "Singleton",]
    cluster_stat = cluster_stat[rev(order(cluster_stat$Freq)),]
    cn = dim(cluster_stat)[1]
    new_membership = paste("Cluster_", 1:cn, sep="")
    cluster_stat = cbind(cluster_stat, new_membership)
    for(ii in 1:cn)
    {  
      levels(c_membership)[levels(c_membership) %in% cluster_stat[ii,1] ] = as.character(cluster_stat[ii,3])
    }
    cluster_stat = as.data.frame(table(c_membership))
    cluster_ID   = as.data.frame(cbind(rownames(cell.cell.jaccard.distance), as.character(c_membership)))
    colnames(cluster_ID) = c("Cell_ID", "GiniClust_membership")
    table(cluster_ID$GiniClust_membership)
    table(c_membership)
    
    #if a cluster smaller than 5% of the total cell number, we call it as a rare cell types cluster.
    cell.num = dim(ExprM.RawCounts)[2] * rare_p    
    rare.cells.list.all = list()
    for(c in levels(c_membership))
    {
      list = rownames(cell.cell.jaccard.distance)[which(c_membership == c)]; 
      if(length(list) < cell.num & c!="Singleton"){
        rare.cells.list.all[[c]]     =  list;
      }
    }
    names(rare.cells.list.all)
    
    clustering_membership_r = data.frame(cbind(rownames(cell.cell.jaccard.distance), as.character(c_membership)))
    colnames(clustering_membership_r) = c("cell.ID","cluster.ID")
    clustering_membership_r[,1] = as.character(clustering_membership_r[,1])
    clustering_membership_r[,2] = as.character(clustering_membership_r[,2])
    write.csv(clustering_membership_r, file=paste(out.folder, "/", exprimentID,"_clusterID.csv", sep=""))
    return(list('cell_cell_dist' = cell.cell.jaccard.distance,'c_membership' = c_membership,
                'clustering_membership_r' = clustering_membership_r,'rare.cell' = rare.cells.list.all))
  }
  
  
  if(data.type == 'qPCR'){
    # correlateion based distance matrix
    ExprM.rgcc.final.log2 <- log2(ExprM.RawCounts.filter[GeneList.final,]+0.1)
    dist.cor <- as.matrix(as.dist(1-cor(ExprM.rgcc.final.log2)))
    cell.cell.cor.distance = as.data.frame(dist.cor+0.01)
    
    title             = paste("eps",eps,sep=".")
    data.mclust       = fpc::dbscan(cell.cell.cor.distance, eps = eps, method = "dist", showplot = FALSE) #MinPts default is 5
    
    o_membership      = factor(paste("db_",data.mclust$cluster ,sep=""))
    if(levels(o_membership)[1]=="db_0"){levels(o_membership)[1] = "Singleton"}
    c_membership = o_membership
    cluster_stat = as.data.frame(table(o_membership))
    cluster_stat = cluster_stat[cluster_stat$o_membership != "Singleton",]
    cluster_stat = cluster_stat[rev(order(cluster_stat$Freq)),]
    cn = dim(cluster_stat)[1]
    new_membership = paste("Cluster_", 1:cn, sep="")
    cluster_stat = cbind(cluster_stat, new_membership)
    for(ii in 1:cn)
    {  
      levels(c_membership)[levels(c_membership) %in% cluster_stat[ii,1] ] = as.character(cluster_stat[ii,3])
    }
    cluster_stat = as.data.frame(table(c_membership))
    cluster_ID   = as.data.frame(cbind(rownames(cell.cell.cor.distance), as.character(c_membership)))
    colnames(cluster_ID) = c("Cell_ID", "GiniClust_membership")
    table(cluster_ID$GiniClust_membership)
    table(c_membership)
    
    #if a cluster smaller than 5% of the total cell number, we call it a rare cell types cluster.
    cell.num = dim(ExprM.RawCounts)[2] * rare_p    #cell.num #[1] 125.45
    rare.cells.list.all = list()
    for(c in levels(c_membership))
    {
      list = rownames(cell.cell.cor.distance)[which(c_membership == c)]; 
      if(length(list) < cell.num & c!="singleton"){
        rare.cells.list.all[[c]]     =  list;
      }
    }
    names(rare.cells.list.all)
    
    clustering_membership_r = data.frame(cbind(rownames(cell.cell.cor.distance), as.character(c_membership)))
    colnames(clustering_membership_r) = c("cell.ID","cluster.ID")
    clustering_membership_r[,1] = as.character(clustering_membership_r[,1])
    clustering_membership_r[,2] = as.character(clustering_membership_r[,2])
    write.csv(clustering_membership_r, file=paste(out.folder, "/", exprimentID,"_clusterID.csv", sep=""))
    return(list('cell_cell_dist' = cell.cell.cor.distance,'c_membership' = c_membership,
                'clustering_membership_r' = clustering_membership_r,'rare.cell' = rare.cells.list.all))
  }
}
####################      End     ####################
