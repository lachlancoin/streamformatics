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

#NEW DB
##OPTIONAL DELETE DB
#.GET("admin/delete", org="Coin", project="RAPIDS1",db="Discovery" ,query=list(delete=T))
.GET("admin/make", org="Coin", project="RAPIDS3",db="Discovery")



##STREAM A FILE
flags = list(format='kallisto')
.POST("dist/stream",org="Coin",project="RAPIDS3",db="Discovery" ,body=list(flags=toJSON(flags)),
      files=list("/home/unimelb.edu.au/lcoin/Data/sAPI/Coin/LCAH/6202_12_S49_R1_001.tsv.gz"))


flags = list(slug_phen="function(x)x", 
             slug_sample = "function(x)toupper(x)",
             id="Sequencing_Sample_ID",
             display=F,
             sep=","
)
.POST("dist/upload_pheno", org="Coin",project="RAPIDS3",db="Discovery",
      body=list(flags=toJSON(flags)),
      files = "/home/unimelb.edu.au/lcoin/Data/sAPI/Coin/LCAH/Validation_meta.csv")

.GET("dist/samples", org="Coin", project="RAPIDS3",db="Discovery",outp="json")
types=.GET("dist/types", org="Coin", project="RAPIDS1",db="Discovery",outp="json")
.GET("dist/variables", org="Coin", project="RAPIDS1",db="Discovery",query = list(type=types[[1]]))


#### FOR DEBUGGING.. COMMANDS EQUIV

setwd("~/github/streamformatics/R")
source("libs.R")
opts = setOptions("~/.sfx")
dbDir="/home/unimelb.edu.au/lcoin/Data/sAPI"

keys = keyEnv$new(dbDir)
user=opts$USER
keys$makeDB(user,"Coin","RAPIDS4","Discovery")

endpoint="dist.R"
all = allEnv$new(dbDir,keys, endpoint)
dist=all$get("Coin","RAPIDS4","Discovery")
user=opts$USER
f="/home/unimelb.edu.au/lcoin/Data/sAPI/Coin/LCAH/6202_12_S49_R1_001.tsv.gz"
sampleID = rev(strsplit(f,"/")[[1]])[1]
flags = list(format="kallisto")
dist$stream(f, sampleID,user, flags)
pheno_f = "/home/unimelb.edu.au/lcoin/Data/sAPI/Coin/LCAH/Discovery_meta.csv"
nmef= rev(strsplit(pheno_f,"/")[[1]])[1]
flags = list(slug_phen="function(x)x", 
             slug_sample = "function(x)toupper(x)",
             id="Sequencing_Sample_ID",
             sep=","
             )
dist$upload_pheno(pheno_f, nmef, user, flags)


