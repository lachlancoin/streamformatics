### SCRIPT FOR TB QUERYING WHO CATALOG

.libPaths("~/R/x86_64-pc-linux-gnu-library/4.1/")
library(curl);library(httr); library(jsonlite);library(readr);library(ggplot2)
source("~/github/streamformatics/R/helper.R")

## move to data directory
setwd("/home/unimelb.edu.au/lcoin/Data/sparsely")


##SET OPTIONS
opts = setOptions("~/.sfx")
if(is.null(opts$USER) || is.null(opts$PASS) || is.null(opts$URL)) stop("!!")
if(opts$URL=="https://api.localhost") httr::set_config(httr::config(ssl_verifypeer = 0L, ssl_verifyhost = 0L))

##TEST API RESPONSIVE
.GET("admin/echo", query = list(msg = "A simple text string"), outp="json", encode="json")

##REGISTER USER to a specific org
.GET("admin/register",org="Coin",query = list(key=opts$API_KEY))


##GET WHO RESISTANCE SNPS
.GET("who/snps",org="Coin", project="WHO", db="2023", outp="tsv")

##RESISTANCE FROM VCF
flags = fromJSON('{"include":["1) Assoc w R","2) Assoc w R - Interim"]}') ## specifies which resistance to show
.POST("who/resistance/vcf",body = list( flags = toJSON(flags)), 
      files=list("ERR10225529.targets.vcf.gz"),
      query=list(sampleID="ERR10225529"),
      org = "Coin",project="WHO",db="2023", outp="tsv", encode="multipart")

.POST("who/resistance/sample", org="Coin", project="WHO",db="2023", query=list(sampleID="ERR10225529"),body = list(flags=toJSON(flags)),outp="tsv")



##CLEAR RESISTANCE FOR SAMPLE
.GET("who/clear",query = list(sampleID="tb-s1.bam"), org="Coin",project="WHO",db="2023", encode="json")

#RESISTANCE FROM BAM
.POST("who/resistance/bam",files=list("tb-s1.bam"),query=list(sampleID="tb-s1.bam"),
      body=list(flags=toJSON(flags)),org="Coin",project="WHO",db="2023",outp="tsv")

##STORED RESISTANCE FOR SAMPLE UNTIL CLEAR
.POST("who/resistance/sample", org="Coin", project="WHO",db="2023", query=list(sampleID="tb-s1.bam"),body = list(flags=toJSON(flags)),outp="tsv")

.POST("who/resistance/sample", org="Coin", project="WHO",db="2023", query=list(sampleID="tb-s1.bam"),body = list(flags=toJSON(flags)),outp="tsv")





###FOLLOWING FOR DEBUGGING PURPOSES - CALLS CLASSES DIRECTLY 
setwd("~/github/streamformatics/R")
source("libs.R")
opts = setOptions("~/.sfx")
dbDir="/home/unimelb.edu.au/lcoin/Data/sAPI"
keys = keyEnv$new(dbDir)
endpoint="who.R"
all = allEnv$new(dbDir,keys, endpoint)
who=all$get("Coin","WHO",db="2023")
user=opts$USER
setwd("/home/unimelb.edu.au/lcoin/Data/sparsely")
flags = list(show=T)


vcf="ERR10225529.vcf.gz"
sampleID=vcf
who$clear(user,sampleID)
aa1=who$vcf(vcf,sampleID, user, flags)
who$resistance(sampleID, user,flags) 


bamf="tb-s1.bam"
bamf="ERR10225529.bam" 
sampleID=bamf
who$clear(user,sampleID)
who$bam(bamf,sampleID,user,flags)
who$resistance(sampleID, user,flags) 




