path=("~/github/sAPI/plumber-thing/")
.libPaths("~/R/x86_64-pc-linux-gnu-library/4.1/")
options(bigmemory.allow.dimnames=TRUE)
library(jsonlite)
library(readr);
library(DBI); 
library(RSQLite);
library(vcfR)
library(plumber)
library(R6)
library(Matrix)
library(tidyr)
library(GenomicRanges)
library(dplyr)
library(readxl)
library(base64enc)
library(Rsamtools); library(GenomicAlignments);
#library(wCorr)
library(confintr)

library(arrow)
library(binom); 
library(pROC); 
library(glmnet)
library(MASS);
library(ggplot2)
library(bigmemory)
library(bigalgebra)
library(biganalytics)
library(data.table)
library(R.utils)

#maxsize=1.0 * 1e9
pathR = paste(path,"R",sep="/")
src1=grep(".R$",dir(pathR,rec=T),v=T)
invisible(try(lapply(paste(pathR,src1,sep="/"), function(x) {print(x);source(x)})))
source("~/github/streamformatics/R/helper.R")
