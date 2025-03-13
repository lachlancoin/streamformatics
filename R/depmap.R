### SCRIPT FOR TRANSCRIPTOMICS

.libPaths("~/R/x86_64-pc-linux-gnu-library/4.1/")
library(curl);library(httr); library(jsonlite);library(readr);library(ggplot2)
source("~/github/streamformatics/R/helper.R")

##SET OPTIONS
opts = setOptions("~/.sfx")
if(is.null(opts$USER) || is.null(opts$PASS) || is.null(opts$URL)) stop("!!")
if(opts$URL=="https://api.localhost") httr::set_config(httr::config(ssl_verifypeer = 0L, ssl_verifyhost = 0L))


##TEST API RESPONSIVE
.GET("admin/echo", query = list(msg = "A simple text string"), outp="json", encode="json")

##REGISTER USER to a specific org
.GET("admin/register",org="Coin",query = list(key=opts$API_KEY))
.GET("admin/list",org=org)


.GET("admin/make", org="Coin", project="Depmap2",db="primary")



##FUNCTION TO PARSE DEPMAP PHENOTYPES
.getDepmapPhenotypes<-function(data_directory, inpfile_, sampfile, treat_file){
  .convert=function(x) gsub("\\.","-",strsplit(x,"\\.\\.")[[1]][1])
  inpfile = paste(data_directory,unlist(inpfile_) ,sep="/")
  names(inpfile) = names(inpfile_)
  sample_file = paste(data_directory,sampfile,sep="/")
  treatment_file = paste(data_directory,treat_file,sep="/")
  samples = read.csv(sample_file,head=T, row.names=1)  
  treatment=read.csv(treatment_file,head=T,row.names=1)
  phen=read.csv(inpfile[['pheno']] ,head=T, row.names=1)
  nme2Treatment=unlist(lapply(dimnames(phen)[[2]],.convert ))
  na_inds =is.na(treatment$name) & !is.na(treatment$broad_id)
  treatment$name[na_inds]=treatment$broad_id[na_inds]
  names(nme2Treatment) =treatment$name[ match(nme2Treatment,treatment$broad_id)]
  incl_samps = dimnames(phen)[[1]]  %in% dimnames(samples)[[1]]
  phen = phen[incl_samps,]
  incl_samps1 = dimnames(samples)[[1]] %in% dimnames(phen)[[1]]
  samples = samples[incl_samps1,]
  ##samples in same order as phens
  samples = samples[match(dimnames(phen)[[1]], dimnames(samples)[[1]]),]
  dimnames(phen)[[2]] = names(nme2Treatment)
  could_not_match = list()
  na_inds=which(is.na(names(nme2Treatment)))
  sampleID = dimnames(phen)[[1]]
  list(drugs = cbind(sampleID, phen), mat2 = cbind(sampleID, samples))
}


dbDir3="~/Data/Depmap"
files = paste(dbDir3,c("drugs_in.csv","phenos_in.csv"),sep="/")
names(files) = files
if(!file.exists(files[[1]])){
  inpfile_=list(pheno="primary-screen-replicate-collapsed-logfold-change.csv")
  mats = .getDepmapPhenotypes(dbDir3, inpfile_, "sample_info1.csv", "primary-screen-replicate-treatment-info.csv")
  write.table(mats$drugs, paste(dbDir,"drugs_in.csv",sep="/"),col.names=T,row.names=F,quote=F,  sep='\t')
  write.table(mats$mat2, paste(dbDir,"phenos_in.csv",sep="/"),col.names=T,row.names=F,quote=F,  sep='\t')
}

flags=list()
flags$sep="\t"
.POST("dist/upload_pheno",org="Coin",project="Depmap2",db="primary" ,body=list(flags=toJSON(flags)),
      files=files)
phens=.GET("dist/phenos", org="Coin", project="Depmap2",db="primary",outp="json")


f= lapply(list( CN="CCLE_gene_cn.csv.gz", expression="CCLE_expression.csv.gz"),function(a) paste(dbDir3,a,sep="/"))
flags$sep=","
flags$type=names(f)
.POST("dist/upload_matrix",org="Coin",project="Depmap2",db="primary" ,query=list(flags=toJSON(flags)),
      files=f)

dist$importMat()
