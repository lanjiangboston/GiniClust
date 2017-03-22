# Author: Lan Jiang
# Contact information: lan_jiang@hms.harvard.edu

#################### Filter data ####################
GiniClust_Filtering <- function(ExprM.RawCounts,out.folder,exprimentID){
  #expressed_cutoff:value lower than this value could be just noise.
  ExpressedinCell_per_gene=apply(ExprM.RawCounts,1,function(x) length(x[x > expressed_cutoff ]))
  nonMir = grep("MIR|Mir", rownames(ExprM.RawCounts), invert = T)  # because Mir gene is usually not accurate 
  Genelist = intersect(rownames(ExprM.RawCounts)[nonMir],rownames(ExprM.RawCounts)[ExpressedinCell_per_gene >= minCellNum])
  ExpressedGene_per_cell=apply(ExprM.RawCounts[Genelist,],2,function(x) length(x[x>0]))
  length(Genelist)
  ExprM.RawCounts.filter = ExprM.RawCounts[Genelist,ExpressedGene_per_cell >= minGeneNum]
  dim(ExprM.RawCounts.filter) 
  write.table(ExprM.RawCounts.filter, file=paste(out.folder, "/", exprimentID, "_gene.expression.matrix.RawCounts.filtered.csv", sep=""), sep=",", row.names = TRUE, col.names = TRUE)
  return(list('raw' = ExprM.RawCounts.filter))
}
####################      End     ####################
