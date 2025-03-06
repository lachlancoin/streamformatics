### SCRIPT FOR running sparsely

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


##PRELOADING THE DATA MATRICES - THIS HAPPENS AUTOMATICALLY, BUT CAN TAKE TIME
.GET("sparsely/load",org="Coin",project="TB",outp="json")

vcf_file = "~/Data/sparsely/biosamples/results/SAMEA111366156/culture/ERR10225541/clair3/ERR10225541.vcf.gz"
cov_file="~/Data/sparsely/biosamples/results/SAMEA111366156/culture/ERR10225541/bam/ERR10225541.bam.cov.gz"
.POST("sparsely/vcf",body = list(),   query = list(sampleID="ERR10225541"), files=list(vcf_file, cov_file), org = "Coin",project="TB", outp="tsv")



flags =list()
flags[['sampleID']] = "tmp1.bam"
#flags$show=T
.POST("sparsely/bam",files=list("~/Data/sparsely/tmp1.bam"),query = list(sampleID="tmp1.bam"),org="Coin",project="TB",outp="tsv")
#flags$sampleID="VIDRL-15139673.bam"
#.POST("sparsely/bam",files=list("~/Data/sparsely/VIDRL-15139673.bam"),body=list(flags=toJSON(flags)),org="Coin",project="TB",outp="tsv")



