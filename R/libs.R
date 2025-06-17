path=paste0(github,"/sAPI/plumber-thing/")
path2=paste0(github,"/fspls2/R")

options(bigmemory.allow.dimnames=TRUE)
library(jsonlite)
library(readr);
library(DBI); 
library(RSQLite);
library(vcfR)
#library(plumber)
library(R6)
library(Matrix)
library(tidyr)
#library(GenomicRanges)
library(dplyr)
#library(readxl)
#library(base64enc)
#library(Rsamtools); library(GenomicAlignments);
#library(wCorr)
library(confintr)

#library(arrow)
#library(binom); 
#library(pROC); 
library(glmnet)
library(MASS);
library(ggplot2)
library(bigmemory)
library(bigalgebra)
library(biganalytics)
library(data.table)
library(ape)
#library(R.utils)

#maxsize=1.0 * 1e9
pathR = paste(path,"R",sep="/")
src1=grep(".R$",dir(pathR,rec=T),v=T)
invisible(try(lapply(paste(pathR,src1,sep="/"), function(x) {print(x);source(x)})))
src1=grep(".R$",dir(path2,rec=T),v=T)
invisible(try(lapply(paste(path2,src1,sep="/"), function(x) {print(x);source(x)})))

source(paste0(github,"/streamformatics/R/helper.R"))
