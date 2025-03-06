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


##OPTIONAL DELETE DB
#.GET("admin/delete", org="Coin", project="RAPIDS1",db="Discovery" ,query=list(delete=T))
.GET("admin/make", org="Coin", project="TB2",db="mtb1")



##STREAM A VCF


.POST("dist/stream",org="Coin",project="TB2",db="mtb1" ,query=list(format="vcf"),
      files=list("~/Data/sparsely/ERR10225529.targets.vcf.gz"))


.GET("dist/samples", org="Coin", project="TB2",db="mtb1",outp="json")
types=.GET("dist/types", org="Coin", project="TB2",db="mtb1",outp="json")
.GET("dist/variables", org="Coin", project="TB2",db="mtb1",query = list(type=types[[1]]))
