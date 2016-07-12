#################### Preprocessing ####################
GiniClust_Preprocess <- function(data.file,data.type,out.folder,exprimentID){
  subdir <- c(out.folder,paste(out.folder,'figures',sep = '/'))
  for(s in subdir){
    if(!file.exists(s)){
      dir.create(s)
    }
  }
  if(data.type == 'RNA-seq'){
    ExprM.RawCounts  <- read.csv(data.file, sep=",", head=T)
    colsum <- apply(ExprM.RawCounts,2,sum)
    ExprM.normCounts <- t(t(ExprM.RawCounts)*CountsForNormalized/colsum)
    write.table(ExprM.RawCounts, file=paste(out.folder, "/", exprimentID, "_rawCounts.csv", sep=""), sep=",", row.names = TRUE, col.names = TRUE, quote=FALSE)
    write.table(ExprM.normCounts, file=paste(out.folder, "/", exprimentID, "_normCounts.csv", sep=""), sep=",", row.names = TRUE, col.names = TRUE, quote=FALSE)
    return(list('raw' = ExprM.RawCounts,'norm' = ExprM.normCounts))
  }
  if(data.type == 'qPCR'){
    ExprM.log2<- read.csv(data.file, sep=",", head=T)
    ExprM.Nor = as.data.frame(2 ^ ExprM.log2 - 1); ##since ExprM is in log2 scale, need to transform back to normal scale, and -1 is important 
    ExprM.RawCounts  = ExprM.Nor
    ExprM.normCounts = ExprM.Nor 
    write.table(ExprM.RawCounts, file=paste(out.folder, "/", exprimentID, "_rawCounts.csv", sep=""), sep=",", row.names = TRUE, col.names = TRUE, quote=FALSE)
    write.table(ExprM.normCounts, file=paste(out.folder, "/", exprimentID, "_normCounts.csv", sep=""), sep=",", row.names = TRUE, col.names = TRUE, quote=FALSE)
    return(list('raw' = ExprM.RawCounts,'norm' = ExprM.normCounts))
  }
}
####################      End     ####################
