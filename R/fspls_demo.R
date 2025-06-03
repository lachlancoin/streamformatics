### SCRIPT FOR TB QUERYING WHO CATALOG

.libPaths("~/R/x86_64-pc-linux-gnu-library/4.1/")
options(bigmemory.allow.dimnames=TRUE)
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
setwd("~/github/streamformatics/R")
source("libs.R")
opts = setOptions("~/.sfx")
dbDir="/home/unimelb.edu.au/lcoin/Data/sAPI"
keys = keyEnv$new(dbDir)
all = allEnv$new(dbDir,keys, "fspls/data.R")

### MAKE NEW DATABAES
all_dist = allEnv$new(dbDir,keys,"dist.R")
keys$makeDB("lcoin","Coin","LCAH","Discovery")
dist1 = all_dist$get("Coin","LCAH","Discovery")
dist2 = all_dist$get("Coin","LCAH","Validation")
flags = list(nrep=5, bigmatrix=F)
datas=all$get1("Coin","LCAH",flags=flags)
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


flags = list(nrep=1,batch=0)
#datas=all$get1("Coin","golub","golub_data",flags=flags)
datas=all$get1("Coin","LCAH",flags=flags, reload=T)

datas$dims()
user=opts$USER

phens = datas$pheno()[1]
phens = "disease_class"
flags = list(pthresh = 5e-5,max=10, topn=20,nrep=5, train="Discovery") ## return can be model, vars or eval
vars = datas$select ( phens, flags)


#vars1 = fromJSON(toJSON(vars))

#vars = read_json("~/Data/sparsely/vars.json")
all_models =datas$makeAllModels(vars,phens,flags)

eval = datas$evaluateAllModels(all_models,phens[1],flags)
.plotEval1(eval,sep="cohort_measure", grid="data~cv", rename=T,shape_color="subpheno", linetype="fullmodel")

predictions =datas$extractPredictions(all_models,phens[1], flags, CV = F);
#aa=roc(predictions[[2]]$y, predictions[[2]]$X0)
.plotArea(predictions, rename=T)


datas$angles(vars,phens,flags)


