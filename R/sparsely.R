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
.GET("admin/list",org="Coin")




##PRELOADING THE DATA MATRICES - THIS HAPPENS AUTOMATICALLY, BUT CAN TAKE TIME SO CAN DO FIRST
flags = list(bigmatrix=FALSE)
.POST("sparsely/load",org="Coin",project="TB",outp="json", body = list(flags = toJSON(flags)))

vcf_file = "~/Data/sparsely/biosamples/results/SAMEA111366156/culture/ERR10225541/clair3/ERR10225541.vcf.gz"
cov_file="~/Data/sparsely/biosamples/results/SAMEA111366156/culture/ERR10225541/bam/ERR10225541.bam.cov.gz"
sampleID="ERR10225541"
.POST("sparsely/vcf",body = list(),   query = list(sampleID=sampleID), files=list(vcf_file, cov_file), org = "Coin",project="TB", outp="tsv")

flags =list()
flags$type="geno"
flags$tool="mykrobe"
.POST("sparsely/sample", org="Coin", project="TB", query=list(sampleID=sampleID),body = list(flags=toJSON(flags)),outp="tsv")


flags[['sampleID']] = "tmp1.bam"
#flags$show=T
.POST("sparsely/bam",files=list("~/Data/sparsely/tmp1.bam"),query = list(sampleID="tmp1.bam"),org="Coin",project="TB",outp="tsv")




#.POST("sparsely/bam",files=list("~/Data/sparsely/ERR10225529.bam"),query = list(sampleID="tmp1.bam"),org="Coin",project="TB",outp="tsv")

#flags$sampleID="VIDRL-15139673.bam"
#.POST("sparsely/bam",files=list("~/Data/sparsely/VIDRL-15139673.bam"),body=list(flags=toJSON(flags)),org="Coin",project="TB",outp="tsv")




###FOLLOWING FOR RUNNING WITHIN R  - THIS IS MORE FOR DEBUGGING
setwd("~/github/streamformatics/R")
source("libs.R")
setwd("~/Data/sparsely")
opts = setOptions("~/.sfx")
dbDir="/home/unimelb.edu.au/lcoin/Data/sAPI"

keys = keyEnv$new(dbDir)
endpoint="sparsely.R"
all = allEnv$new(dbDir,keys, endpoint)
msim=all$get1("Coin","TB",useBigMatrix = F)
user=opts$USER

flags = list(type="pheno", tool="mykrobe")

vcf_file = "ERR10225529.vcf.gz"
cov_file="ERR10225529.bam.cov.gz" 
bamf="ERR10225529.bam"


cov_file="barcode03.cov.gz"
vcf_file="barcode03.targets.vcf.gz"
sampleID=vcf_file
msim$vcf(vcf_file, cov_file, sampleID, user)
flags$type="pheno"
msim$get_resistance(sampleID, user,flags)
#msim$get_resistance(sampleID, user, flags)

bamf="ERR10225529.bam"
bamf="barcode03.bam"
bamf="tb-s1.bam"
bamf="tmp1.bam"
sampleID=bamf
flags$max_reads=50000
flags$step=5000
flags$show_step=T
msim$bam(bamf, sampleID, user, flags)
flags$type="pheno"
msim$get_resistance(sampleID, user, flags)

