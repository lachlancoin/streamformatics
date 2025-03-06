### SCRIPT FOR TRANSCRIPTOMICS

.libPaths("~/R/x86_64-pc-linux-gnu-library/4.1/")
library(curl);library(httr); library(jsonlite);library(readr);library(ggplot2)
source("~/github/streamformatics/R/helper.R")

##SET OPTIONS
opts = setOptions("~/.sfx")
if(is.null(opts$USER) || is.null(opts$PASS) || is.null(opts$URL)) stop("!!")
if(opts$URL=="https://api.localhost") httr::set_config(httr::config(ssl_verifypeer = 0L, ssl_verifyhost = 0L))

##TEST API RESPONSIVE
a=.GET("admin/echo", query = list(msg = "A simple text string"), outp="json", encode="json")
print(a)

##REGISTER USER to a specific org
.GET("admin/register",org="Coin",query = list(key=opts$API_KEY))
.GET("admin/list",org=org)

#NEW DB
##OPTIONAL DELETE DB
#.GET("admin/delete", org="Coin", project="RAPIDS1",db="Discovery" ,query=list(delete=T))
.GET("admin/make", org="Coin", project="RAPIDS1",db="Discovery")



##STREAM A FILE
.POST("dist/stream",org="Coin",project="RAPIDS1",db="Discovery" ,query=list(format="kallisto"),
      files=list("/home/unimelb.edu.au/lcoin/Data/sAPI/Coin/LCAH/6202_12_S49_R1_001.tsv.gz"))

.GET("dist/samples", org="Coin", project="RAPIDS1",db="Discovery",outp="json")
types=.GET("dist/types", org="Coin", project="RAPIDS1",db="Discovery",outp="json")
.GET("dist/variables", org="Coin", project="RAPIDS1",db="Discovery",query = list(type=types[[1]]))


