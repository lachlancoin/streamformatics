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


##GET WHO RESISTANCE SNPS
.GET("who/snps",org="Coin", project="WHO", db="2023", outp="tsv")

##RESISTANCE FROM VCF
flags = fromJSON('{"include":["1) Assoc w R","2) Assoc w R - Interim"]}')
.POST("who/resistance/vcf",body = list( flags = toJSON(flags)), 
      files=list("~/Data/sparsely/ERR10225529.targets.vcf.gz"),
      query=list(sampleID="ERR10225529"),
      org = "Coin",project="WHO",db="2023", outp="tsv", encode="multipart")

#RESISTANCE FROM BAM
flags$show=TRUE
.POST("who/resistance/bam",files=list("~/Data/sparsely/tmp1.bam"),
      query=list(sampleID="tmp1.bam"),
      body=list(flags=toJSON(flags)),org="Coin",project="WHO",db="2023",outp="tsv")

##STORED RESISTANCE FOR SAMPLE UNTIL CLEAR
.POST("who/resistance/sample", org="Coin", project="WHO",db="2023", query=list(sampleID="tmp1.bam"),body = list(flags=toJSON(flags)),outp="tsv")

#CLEAR CACHE FOR A SAMPLE
.GET("who/clear",query = list(sampleID="tmp1.bam"), org="Coin",project="WHO",db="2023", encode="json")


