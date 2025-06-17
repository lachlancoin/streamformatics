### SCRIPT FOR TB QUERYING WHO CATALOG
rlib="~/R/x86_64-pc-linux-gnu-library/4.4/"
.libPaths(rlib)

github="/data/gpfs/projects/punim1068/github"
options(bigmemory.allow.dimnames=TRUE)
#install.packages("readr",lib=rlib)
library(curl);library(httr); library(jsonlite);library(ggplot2)
library(readr);
source(paste0(github,"/streamformatics/R/helper.R"))

##SET OPTIONS
opts = setOptions("~/.sfx")
if(is.null(opts$USER) || is.null(opts$PASS) || is.null(opts$URL)) stop("!!")
if(opts$URL=="https://api.localhost") httr::set_config(httr::config(ssl_verifypeer = 0L, ssl_verifyhost = 0L))

##TEST API RESPONSIVE
.GET("admin/echo", query = list(msg = "A simple text string"), outp="json", encode="json")

##REGISTER USER to a specific org
.GET("admin/register",org="Coin",query = list(key=opts$API_KEY))
.GET("admin/make", org="Coin", project="golub",db="golub_data")


#DOWNLOAD RDS FILE FROM GITHUB
url="https://github.com/dn-ra/FSPLS-publication-repo/blob/master/input/golub_data/golub.prepd.Rds"
download_file = "./golub_prepd.Rds"
GET(url, write_disk(download_file, overwrite=TRUE))

#inputs = as.list(grep(".Rds",dir("~/github/FSPLS-publication-repo/input",rec=T,full=T),v=T))
#names(inputs)=lapply(inputs, function(str)rev(strsplit(str,"/")[[1]])[2])
db="golub_data"
#.GET("admin/delete", org="Coin", project="golub", db="golub",query=list(delete=T))
files = list("golub_data" = download_file)
.POST("dist/upload_rds", org="Coin", project="golub",db=db,files = inputs[db], body = list(flags=toJSON(list(force=T))),outp="json")


flags = list(bigmatrix=FALSE,reload=T)
.POST("fspls/data/load",org="Coin",project="golub",outp="json", body = list(flags=toJSON(flags)))


phens=.GET("fspls/data/pheno", org="Coin", project="golub",outp="json")

flags = list(pthresh = 1e-2, topn=10,beam=1,train="golub_data",test="golub_data") ## return can be model, vars or eval
vars=.POST("fspls/data/train",org="Coin",project="golub", body= list(flags = toJSON(flags), phens =toJSON(phens)), outp="json")
print(vars)
models=.POST("fspls/data/makeModels",org="Coin",project="golub", body= list(flags = toJSON(flags), phens =toJSON(phens),vars = toJSON(vars)), outp="json")


eval =.POST("fspls/data/evaluate",org="Coin",project="golub", body= list(flags = toJSON(flags), phens =toJSON(phens),models = toJSON(models)), outp="tsv")

print(eval)

vars[[2]] = vars[[1]][1]
names(vars) = lapply(vars, function(v) paste(names(v), collapse="."))
angles=.POST("fspls/data/angles",org="Coin",project="golub", outp="json",  body= list(flags = toJSON(flags), phens =toJSON(phens), vars = toJSON(vars)))
print(angles)

#vars[[2]] = vars[[1]][length(vars[[1]])]
#vars[[3]] = vars[[1]][length(vars[[1]])-1]
#vars[[4]] = vars[[1]][1:(length(vars[[1]])-1)]
names(vars) = lapply(vars, function(v) paste(names(v), collapse="."))
pvalues = .POST("fspls/data/pvalue",org="Coin",project="golub", outp="json",  body= list(flags = toJSON(flags), phens =toJSON(phens), vars = toJSON(vars)))
print(pvalues)



ggplot(eval)+geom_point(aes(x=numvars, y=-value, shape=fullmodel, color=data))+
  geom_line(aes(x=numvars, y=-value, linetype=fullmodel, color=data))+
  facet_grid("pheno~type")







###DEBUGGING

setwd(paste0(github,"/streamformatics/R"))
source("libs.R")
opts = setOptions("~/.sfx")
dbDir="/home/lcoin/punim1140/Depmap_LCOIN/sAPI"
keys = keyEnv$new(dbDir)
all = allEnv$new(dbDir,keys, "fspls/data.R")

### MAKE NEW DATABAES
all_dist = allEnv$new(dbDir,keys,"dist.R")
keys$makeDB("lcoin","Coin","LCAH","Discovery")

###DEPMAP
#keys$dbs("Coin","Depmap")
#keys$removeDB(opts$USER, "Coin","Depmap","primary",deleteDir=T)
#keys$makeDB(opts$USER, "Coin","Depmap","primary")

dist1 = all_dist$get("Coin","Depmap","primary",reset=F)
dist1$sampleID()

##UPLOADING DEPMAP
depmap="/home/unimelb.edu.au/lcoin/Data/Depmap"
filenames=list(rna="CCLE_expression.csv.gz",cn="CCLE_gene_cn.csv.gz"  )
flags = list()
lapply(names(filenames), function(nme){
  print(nme)
    mat = data.frame(data.table::fread(paste(depmap,filenames[[nme]],sep="/"),sep=","))
    rownames(mat) = mat[,1]
    colnames(mat) = unlist(lapply(colnames(mat), function(str) strsplit(str,"\\..")[[1]][1]))
    #dimnames(mat)[[1]] = mat[,1]
    print("importing")
    dist1$importMat(mat[,-1], filenames[[nme]], opts$USER, flags, type=nme)
    print("done")
})
nme1 = "primary-screen-replicate-collapsed-logfold-change.csv.gz"
filenme1 = paste(depmap,nme1,sep="/")
flags = list()
flags[['sep']]=','; flags[['id']]="V2"
phenos1 = dist1$getPheno("drugs")
phenos1$upload_pheno(filenme1,nme1,opts$USER,flags, samples=dist1$sampleID())
phenos2 = dist1$getPheno("cells")
nme2 = "sample_info.csv.gz"
filenme2 = paste(depmap,nme2,sep="/")
flags$slug_sample = '["function(x) x"]'
flags$id = "DepMap_ID";flags$sep=","
phenos2$upload_pheno(filenme2, nme2, opts$USER, flags, samples = dist1$sampleID())
#####

flags = list(nrep=0,batch=1,bigmatrix=F,splitBy="disease",nphenos=50000, keep="pancreatic")
datasAll = all$get1("Coin","Depmap",flags=flags, reload=F)



dist1 = all_dist$get("Coin","LCAH","Discovery")
dist2 = all_dist$get("Coin","LCAH","Validation")
flags = list(nrep=5, bigmatrix=F, nphenos=5)
datasAll=all$get1("Coin","LCAH",flags=flags)
#datas2=all$get1("Coin","LCAH","Validation",flags=flags)
dir="/home/unimelb.edu.au/lcoin/Data/sAPI/Coin/LCAH/LCAH_data"
dir1="/home/unimelb.edu.au/lcoin/Data/sAPI/Coin/LCAH"
nme="22022024_RAPIDS_Pheonix_scores.csv"
nme1 = "Validation_meta.csv"
filenme = paste(dir,nme,sep="/")
filenme1 = paste(dir1,nme1,sep="/")

phenos = dist1$phenos;
phenos2 = dist2$phenos
phenos$upload_pheno(filenme1,nme1,opts$USER,flags, samples=dist1$sampleID())
phenos2$upload_pheno(filenme1,nme1,opts$USER,flags, samples=dist2$sampleID())

phenos$upload_pheno(filenme1,nme1,opts$USER,flags, samples=samples)
export = phenos$exportPheno()
#dist1$upload_pheno()
#####


#flags = list(nrep=10,batch=0)
#datasAll=all$get1("Coin","golub","golub_data",flags=flags)
#datasAll=all$get1("Coin","LCAH",flags=flags, reload=T)

##RUN
datasAll$dims()
user=opts$USER

phens = datasAll$pheno(sep=T)
#grep(bestpheno,unlist(phens))
#phens = phens[4499]
options("fspls.types"=
          fromJSON('{"gaussian": ["correlation","var","mad"],"binomial":"AUC","multinomial":"AUC_all","ordinal" : "AUC_all"}'))
#phens = phens[1:4]
flags1 = list(var_thresh = 0.05, 
              pthresh = 1e-6,max=10, topn=20,nrep=5,beam=1, train=grep("Pancreatic", names(datasAll$datas),v=T)) ## return can be model, vars or eval
flags1$test=flags1$train

  vars = datasAll$select ( phens, flags1,verbose=T)
  if(length(names(vars))==0) return(NULL)

  full_inds= lapply(vars$inds, function(vars1) vars1[unlist(lapply(vars1, function(var) "full" %in% names(var)))])
  full_inds = full_inds[unlist(lapply(full_inds, length))>0]
  print(full_inds)
  
  all_models =datasAll$makeAllModels(vars,phens,flags1)

  eval0 = datasAll$evaluateAllModels(all_models,phens,flags1)
  eval0 = subset(eval0, !is.na(mid))
  
  toSave = list(vars = vars, all_models = all_models, eval0 = eval0)  
  saveRDS(toSave,"all1.rds")
  eval_cv_full = subset(eval0,isfull & cv & measure=="correlation" ) 
  hist(eval_cv_full$mid,br=100)
  
  eval_cv = (subset(eval0,cv==T & measure=="var"))
 # which.max(eval_cv_full$mid)
  eval1 = .calcEval1(eval0,sep_by="pheno")
  #eval1 = eval1[grep("avg", eval1$cv, inv=T),]
  eval1 = subset(eval1, cv_full=="CV= TRUE FULL= TRUE")
  ggp=.plotEval1(eval1, grid="measure~cv_full", rename=T,shape_color="pheno", linetype="fullmodel",logy=F, legend=F,showranges=F)
  eval1
#})
eval0 = .merge1_new(eval_all)
eval2 = subset(eval0,cv=="CV= avg" & measure=="var")
#eval2 = subset(eval_,cv=="CV= avg" & measure=="var")
topmid =tail(sort(eval2$mid,decr=T),1)
best_pheno=unique(eval2$pheno[match(topmid,eval2$mid)])

eval3 = subset(eval0,pheno %in% best_pheno)
eval3$pheno=factor(eval3$pheno, levels = best_pheno)
#eval2 = subset(eval,pheno==eval$p)


ggps=.plotEval1(eval3, grid="pheno~cv", rename=T,shape_color="measure", linetype="fullmodel",logy=T)
ggps
#ggps=.plotEval1(eval,sep="pheno", grid="data~cv", rename=T,shape_color="subpheno", linetype="fullmodel")

predictions =datasAll$extractPredictions(all_models,phens[1], flags, CV = F);
#aa=roc(predictions[[2]]$y, predictions[[2]]$X0)
.plotArea(predictions, rename=T)


datasAll$angles(vars,phens,flags)


