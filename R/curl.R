.libPaths("~/R/x86_64-pc-linux-gnu-library/4.1/")
library(curl);library(httr); library(jsonlite);library(readr);library(ggplot2)

setOptions<-function(config){
    if(!file.exists(config))stop("need config file specifying USER,PASS and KEY")
    config = t(read.csv(config,sep="=",header=F))
    config1 = as.list(config[2,])
    names(config1) = config[1,]
    config1 = rev(config1)
    config1 = config1[!duplicated(names(config1))]
    options(config1)
    config1
}
opts = setOptions("~/.sparsely")
if(is.null(opts$USER) || is.null(opts$PASS) || is.null(opts$URL)) stop("!!")
if(opts$URL=="https://api.localhost") httr::set_config(httr::config(ssl_verifypeer = 0L, ssl_verifyhost = 0L))



.GET<-function(endpoint, org=NULL, project=NULL, db=NULL, query = list(),body = list(), outp="json", encode="json"){
  auth = authenticate(getOption("USER"), getOption("PASS"));baseurl = getOption("URL")
  url=paste(baseurl,endpoint,sep="/")
  if(!is.null(org)) url = paste(url, org,sep="/")
  if(!is.null(project)) url = paste(url,project,sep="/")
  if(!is.null(db)) url = paste(url,db,sep="/")
  print(url)
  #body = list()
  #if(length(flags)>0) body = list(flags = toJSON(flags))
  print(body)
  a=GET(url, query = query,body=body,encode=encode,auth)
  a1=rawToChar(a$content)
  if(outp=="json") return(fromJSON(a1))
  if(outp=="tsv") return(readr::read_tsv(a1))
  return(a1)
}
.POST<-function(endpoint, org=NULL, project=NULL, db=NULL, query=list(),body = list(), files=NULL, outp="json",encode="multipart"){
  auth = authenticate(getOption("USER"), getOption("PASS"));baseurl = getOption("URL")
  url=paste(baseurl,endpoint,sep="/")
  if(!is.null(org)) url = paste(url, org,sep="/")
  if(!is.null(project)) url = paste(url,project,sep="/")
  if(!is.null(db)) url = paste(url,db,sep="/")
  print(url)
#  body = list(flags=toJSON(flags))
  if(length(files)>0){
    if(length(files)>1){
    names(files) = rep("upload", length(files))
    uploads=lapply(files, function(f) upload_file(f))
    body = c(body, uploads)
    }else{
      body$upload =upload_file(files[[1]])
    }
  }
  if(length(query)>0){
  a=POST(url,body = body , query=query,auth,encode=encode) #,encode="multipart")
  }else{
    a = POST(url, body = body, query = query, auth, encode=encode)
  }
  a1=rawToChar(a$content)
  if(outp=="json") return(fromJSON(a1))
  if(outp=="tsv") return(readr::read_tsv(a1))
  return(a1)
}



######################  ADMIN ednpoints ##################33
a=.GET("admin/echo", query = list(msg = "A simple text string"), outp="json", encode="json")
print(a)

##REGISTER USER to a specific org
.GET("admin/register",org="Coin",query = list(key=opts$API_KEY))

.GET("admin/list",org=org)
.GET("admin/list", org=org, project="pdig")

.GET("admin/make", org="Coin", project="TB1")
.GET("admin/make", org="Coin", db="mtb1",project="TB1")





###################################################################  WHO ENDPONTS #####

##GET WHO RESISTANCE SNPS
.GET("who/snps",org="Coin", project="WHO", db="2023", outp="tsv")

##RESISTANCE FROM VCF
flags = fromJSON('{"include":["1) Assoc w R","2) Assoc w R - Interim"]}')
flags[['sampleID']]="ERR10225529"
.POST("who/resistance/vcf",body = list( flags = toJSON(flags)), 
      files=list("~/Data/sparsely/ERR10225529.targets.vcf.gz"),
      org = "Coin",project="WHO",db="2023", outp="tsv", encode="multipart")
#RESISTANCE FROM BAM
flags$show=TRUE
flags[['sampleID']] = "tmp1.bam"
.POST("who/resistance/bam",files=list("~/Data/sparsely/tmp1.bam"),body=list(flags=toJSON(flags)),org="Coin",project="WHO",db="2023",outp="tsv")

##STORED RESISTANCE FOR SAMPLE UNTIL CLEAR
.POST("who/resistance/sample", org="Coin", project="WHO",db="2023", body = list(flags=toJSON(flags)),outp="tsv")

#CLEAR CACHE FOR A SAMPLE
.GET("who/clear",query = list(sampleID=flags$sampleID), org="Coin",project="WHO",db="2023", encode="json")




################################################################   DIST ENDPOINTS - FOR UPLOADING DATA
##Make a new database
.GET("admin/make", org="Coin", project="pdig")



inputs = as.list(grep(".Rds",dir("~/github/FSPLS-publication-repo/input",rec=T,full=T),v=T))
names(inputs)=lapply(inputs, function(str)rev(strsplit(str,"/")[[1]])[2])
db="golub_data"
.GET("admin/delete", org="Coin", project="pdig", db=db,query=list(delete=T))
.GET("admin/make", org="Coin", project="pdig", db=db)
#.GET("admin/delete", org="Coin", project="pdig2",query=list(delete=T))
.POST("dist/upload_rds", org="Coin", project="pdig",db=db,files = inputs[db], body = list(flags=toJSON(list(force=T))),outp="json")

.GET("admin/delete", org="Coin", project="TB1", db="mtb1",query=list(delete=T))
.GET("admin/make", org="Coin", project="TB7", db="mtb1")

.POST("dist/stream",org="Coin",project="TB7",db="mtb1" ,query=list(format="vcf"),
      files=list("~/Data/sparsely/ERR10225529.targets.vcf.gz"))


file1="/home/unimelb.edu.au/lcoin/Data/sAPI/Coin/LCAH/6202_12_S49_R1_001.tsv.gz"

.POST("dist/stream",org="Coin",project="RAPIDS",db="Discovery" ,query=list(format="kallisto"),
      files=list(file1))

#CHECK
.GET("dist/types", org="Coin", project="pdig",db=db,outp="json")
.GET("dist/samples", org="Coin", project="pdig",db=db,outp="json")
.GET("dist/phenos", org="Coin", project="pdig",db=db,outp="json")
.GET("dist/cats", org="Coin", project="pdig",db=db,outp="json", query=list(pheno="x"))
.GET("dist/variables", org="Coin", project="pdig",db=db,outp="json",query = list(type="rna"))


#######################################################    SPARSELY END PONTS ############3
flags =list()
vcf_file = "~/Data/sparsely/biosamples/results/SAMEA111366156/culture/ERR10225541/clair3/ERR10225541.vcf.gz"
cov_file="~/Data/sparsely/biosamples/results/SAMEA111366156/culture/ERR10225541/bam/ERR10225541.bam.cov.gz"

##THIS DOESNT WORK VIA R (due to uploading multiple files BUT DOES WORK VIA SSH
#.POST("sparsely/vcf",body = list( flags = toJSON(flags)),  file=list(vcf_file, cov_file), org = "Coin",project="TB", outp="json")
     
     
##PRELOADING THE DATA MATRICES
.GET("sparsely/load",org="Coin",project="TB",outp="json")

flags[['sampleID']] = "tmp1.bam"
#flags$show=T
.POST("sparsely/bam",files=list("~/Data/sparsely/tmp1.bam"),query = list(sampleID="tmp1.bam"),org="Coin",project="TB",outp="tsv")
#flags$sampleID="VIDRL-15139673.bam"
#.POST("sparsely/bam",files=list("~/Data/sparsely/VIDRL-15139673.bam"),body=list(flags=toJSON(flags)),org="Coin",project="TB",outp="tsv")

flags[['sampleID']] = "ERR10225541"

.POST("sparsely/vcf",body = list(),   query = list(sampleID="ERR10225541"), files=list(vcf_file, cov_file), org = "Coin",project="TB", outp="tsv")
    
     



#######################################################    FSPLS END PONTS ############3
.GET("fspls/data/pheno", org="Coin", project="pdig",outp="json")
#.GET("fspls/data/categories", org="Coin", project="pdig",outp="json")



phens = fromJSON('["x.x"]')
flags = list(return="vars",pthresh = 5e-2, topn=10,beam=1) ## return can be model, vars or eval
vars=.GET("fspls/data/train",org="Coin",project="pdig", query= list(flags = toJSON(flags), phens =toJSON(phens),vars = toJSON(vars)), outp="json")
flags$return="model"
models=.GET("fspls/data/train",org="Coin",project="pdig", query= list(flags = toJSON(flags), phens =toJSON(phens),vars = toJSON(vars)), outp="json")


eval =.POST("fspls/data/evaluate",org="Coin",project="pdig", query= list(flags = toJSON(flags), phens =toJSON(phens)),body = list(models = toJSON(models)), outp="tsv")
eval1 =.POST("fspls/data/evaluate",org="Coin",project="pdig", query= list(flags = toJSON(flags), phens =toJSON(phens)),body = list(vars = toJSON(vars)), outp="tsv")

print(eval)
print(eval1)

vars[[2]] = vars[[1]][1]
names(vars) = lapply(vars, function(v) paste(names(v), collapse="."))

angles=.GET("fspls/data/angles",org="Coin",project="pdig", outp="json",  query= list(flags = toJSON(flags), phens =toJSON(phens), vars = toJSON(vars)))


#vars[[2]] = vars[[1]][length(vars[[1]])]
#vars[[3]] = vars[[1]][length(vars[[1]])-1]
#vars[[4]] = vars[[1]][1:(length(vars[[1]])-1)]
names(vars) = lapply(vars, function(v) paste(names(v), collapse="."))
res = .GET("fspls/data/pvalue",org="Coin",project="pdig", outp="json",  query= list(flags = toJSON(flags), phens =toJSON(phens), vars = toJSON(vars)))




ggplot(eval)+geom_point(aes(x=numvars, y=-value, shape=fullmodel, color=data))+
  geom_line(aes(x=numvars, y=-value, linetype=fullmodel, color=data))+
  facet_grid("pheno~type")

      
      