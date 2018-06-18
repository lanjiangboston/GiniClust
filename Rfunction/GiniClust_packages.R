# Author: Lan Jiang
# Contact information: lan_jiang@hms.harvard.edu

#########################################################
### A) Installing and loading required packages
#########################################################

if (!suppressWarnings(require("ggplot2",quietly = TRUE))) {
   install.packages("ggplot2", dependencies = TRUE, repos="https://cran.r-project.org")
   }
if (!suppressWarnings(require("hydroGOF",quietly = TRUE))) {
   install.packages("hydroGOF", dependencies = TRUE, repos="https://cran.r-project.org")
   }
if (!suppressWarnings(require("glmnet",quietly = TRUE))) {
   install.packages("glmnet", dependencies = TRUE, repos="https://cran.r-project.org")
   }
if (!suppressWarnings(require("StatMatch",quietly = TRUE))) {
   install.packages("StatMatch", dependencies = TRUE, repos="https://cran.r-project.org")
   }
if (!suppressWarnings(require("Rtsne",quietly = TRUE))) {
   install.packages("Rtsne", dependencies = TRUE, repos="https://cran.r-project.org")
   }
if (!suppressWarnings(require("fpc",quietly = TRUE))) {
   install.packages("fpc", dependencies = TRUE, repos="https://cran.r-project.org")
   }
if (!suppressWarnings(require("GA",quietly = TRUE))) {
   install.packages("GA", dependencies = TRUE, repos="https://cran.r-project.org")
   }
if (!suppressWarnings(require("MASS",quietly = TRUE))) {
   install.packages("MASS", dependencies = TRUE, repos="https://cran.r-project.org")
   }
if (!suppressWarnings(require("session",quietly = TRUE))) {
   install.packages("session", dependencies = TRUE, repos="https://cran.r-project.org")
   }
if (!suppressWarnings(require("Matrix",quietly = TRUE))) {
   install.packages("Matrix", dependencies = TRUE, repos="https://cran.r-project.org")
   }
if (!suppressWarnings(require("vegan",quietly = TRUE))) {
   install.packages("vegan", dependencies = TRUE, repos="https://cran.r-project.org")
   }
if (!suppressWarnings(require("data.table",quietly = TRUE))) {
   install.packages("data.table", dependencies = TRUE, repos="https://cran.r-project.org")
   }
if (!suppressWarnings(require("reshape",quietly = TRUE))) {
   install.packages("reshape", dependencies = TRUE, repos="https://cran.r-project.org")
   }
if (!suppressWarnings(require("SparseM",quietly = TRUE))) {
   install.packages("SparseM",dependencies = TRUE, repos="https://cran.r-project.org")
   }
if (!suppressWarnings(require("abind",quietly = TRUE))) {
   install.packages("abind", dependencies = TRUE, repos="https://cran.r-project.org")
   }
if (!suppressWarnings(require("drc",quietly = TRUE))) {
   install.packages("drc", dependencies = TRUE, repos="https://cran.r-project.org")
   }
if (!suppressWarnings(require("gplots",quietly = TRUE))) {
   install.packages("gplots",dependencies = TRUE, repos="https://cran.r-project.org")
}


#Load libraries
library(ggplot2)
library(hydroGOF)
library(glmnet)
library(StatMatch)
library(Rtsne)
library(fpc)
library(GA)
library(MASS)
library(session)
library(Matrix)
library(vegan)
library(data.table)
library(reshape)
library(abind)
library(drc)
