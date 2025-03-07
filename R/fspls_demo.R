### SCRIPT FOR TB QUERYING WHO CATALOG

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





phens=.GET("fspls/data/pheno", org="Coin", project="golub",outp="json")

flags = list(return="vars",pthresh = 5e-2, topn=10,beam=1,train="golub_data",test="golub_data") ## return can be model, vars or eval
vars=.GET("fspls/data/train",org="Coin",project="golub", query= list(flags = toJSON(flags), phens =toJSON(phens),vars = toJSON(vars)), outp="json")
print(vars)
flags$return="model"
models=.GET("fspls/data/train",org="Coin",project="golub", query= list(flags = toJSON(flags), phens =toJSON(phens),vars = toJSON(vars)), outp="json")


eval =.POST("fspls/data/evaluate",org="Coin",project="golub", query= list(flags = toJSON(flags), phens =toJSON(phens)),body = list(models = toJSON(models)), outp="tsv")
eval1 =.POST("fspls/data/evaluate",org="Coin",project="golub", query= list(flags = toJSON(flags), phens =toJSON(phens)),body = list(vars = toJSON(vars)), outp="tsv")

print(eval)
print(eval1)

vars[[2]] = vars[[1]][1]
names(vars) = lapply(vars, function(v) paste(names(v), collapse="."))
angles=.GET("fspls/data/angles",org="Coin",project="golub", outp="json",  query= list(flags = toJSON(flags), phens =toJSON(phens), vars = toJSON(vars)))
print(angles)

#vars[[2]] = vars[[1]][length(vars[[1]])]
#vars[[3]] = vars[[1]][length(vars[[1]])-1]
#vars[[4]] = vars[[1]][1:(length(vars[[1]])-1)]
names(vars) = lapply(vars, function(v) paste(names(v), collapse="."))
pvalues = .GET("fspls/data/pvalue",org="Coin",project="golub", outp="json",  query= list(flags = toJSON(flags), phens =toJSON(phens), vars = toJSON(vars)))
print(pvalues)



ggplot(eval)+geom_point(aes(x=numvars, y=-value, shape=fullmodel, color=data))+
  geom_line(aes(x=numvars, y=-value, linetype=fullmodel, color=data))+
  facet_grid("pheno~type")
