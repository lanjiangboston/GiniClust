# Authors: Lan Jiang; Huidong Chen
# Contact information: lan_jiang@hms.harvard.edu; 
#                      hdchen@jimmy.harvard.edu; 

############################### work directory ############################### 
workdir = getwd() 
setwd(workdir)
if(!file.exists('Library')){
  dir.create('Library')
}

############################### Command Interface ############################### 
.libPaths(c(paste(workdir,'/Library',sep=""),.libPaths()))

if (!suppressWarnings(require("optparse",quietly = TRUE))) {
  install.packages("optparse",dependencies = TRUE, repos="http://cran.r-project.org")
}
library(optparse,quietly = TRUE)

option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="input dataset file name", metavar="character"),
  make_option(c("-t", "--type"), type="character", default=NULL, 
              help="input dataset type: choose from 'qPCR' or 'RNA-seq' ", metavar="character"),  
  make_option(c("-o", "--out"), type="character", default="results", 
              help="output folder name [default= %default]", metavar="character"),
  make_option(c("-e", "--epsilon"), type="double", default=NULL, 
            help="DBSCAN epsilon parameter qPCR:[default=0.25],RNA-seq:[default=0.5]", metavar="double"),
  make_option(c("-m", "--minPts"), type="integer", default=NULL, 
              help="DBSCAN minPts parameter qPCR:[default=5],RNA-seq:[default=3]", metavar="integer")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$file)){
  print_help(opt_parser)
  stop("Input dataset file name must be supplied", call.=FALSE)
}else if(!file.exists(opt$file)){
  print_help(opt_parser)
  stop("Input dataset file does not exist", call.=FALSE)
}else if (is.null(opt$type)){
  print_help(opt_parser)
  stop("Input dataset type must be supplied", call.=FALSE)
}else if (!(opt$type %in% c('qPCR','RNA-seq'))){
  print_help(opt_parser)
  stop("Input dataset type must be chosen from 'qPCR' or 'RNA-seq' ", call.=FALSE)
}

data.file = opt$file
data.type = opt$type
out.folder = opt$out
exprimentID = basename(data.file)


############################### GiniClust pipeline ###############################

#parameters and installation of packages
source("Rfunction/GiniClust_parameters.R")
source("Rfunction/GiniClust_packages.R")

#preprocess
source("Rfunction/GiniClust_Preprocess.R")
ExprM.Results = GiniClust_Preprocess(data.file,data.type,out.folder,exprimentID)
ExprM.RawCounts = ExprM.Results$raw

#filtering
source("Rfunction/GiniClust_Filtering.R")
ExprM.Results.filter = GiniClust_Filtering(ExprM.RawCounts,ExprM.normCounts,out.folder,exprimentID)
ExprM.RawCounts.filter = ExprM.Results.filter$raw

#gene selection
source("Rfunction/GiniClust_Fitting.R")
GeneList.final = GiniClust_Fitting(data.type,ExprM.RawCounts.filter,out.folder,exprimentID)

#clustering
source("Rfunction/GiniClust_Clustering.R")
Cluster.Results = GiniClust_Clustering(data.type,ExprM.RawCounts.filter,GeneList.final,eps,MinPts,out.folder,exprimentID)
cell.cell.distance = Cluster.Results$cell_cell_dist
c_membership = Cluster.Results$c_membership
clustering_membership_r = Cluster.Results$clustering_membership_r
rare.cells.list.all = Cluster.Results$rare.cell

#tSNE visualzation
source("Rfunction/GiniClust_tSNE.R")
GiniClust_tSNE(data.type,c_membership,cell.cell.distance,perplexity,out.folder,exprimentID)

#check the clustering results
table(c_membership)
print(rare.cells.list.all)

#Analysis of differentially expressed genes 
if(data.type == 'RNA-seq'){
  source("Rfunction/DE_MAST.R")
  DE_MAST(ExprM.RawCounts.filter,rare.cells.list.all,out.folder,exprimentID)
}
if(data.type == 'qPCR'){
  source("Rfunction/DE_t_test.R")
  DE_t_test(ExprM.RawCounts.filter,rare.cells.list.all,c_membership,out.folder,exprimentID)
}
