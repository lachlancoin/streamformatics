


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

.plotAll<-function(evals_all, max_len =16){
  range = seq(max(1, length(evals_all)-max_len) ,length(evals_all))
  evals2=.merge1_new(evals_all[range])
  evals2$pheno=dEnv$convertFromCode(evals2$pheno)
  ggps=.plotEval1(evals2,legend=T, grid0=c("pheno"), grid1=NULL, #"measure",
                  shape_color=c("cv_full"),
                  sep_by=c("data", "measure"),
                  showranges=T,linetype="cv_full",
                  scales="free",title =names(phens)[1], title1=names(phens1[[1]]))
  ggps
}

opts = setOptions("~/.sfx")
rlib = opts$rlib
dbDir = opts$dbDir
github=opts$github
.libPaths(rlib)

options(bigmemory.allow.dimnames=TRUE)
#install.packages("readr",lib=rlib)
library(curl);library(httr); library(jsonlite);library(ggplot2)
library(readr); library(cowplot); library(ggrepel)

##SET OPTIONS
if(is.null(opts$USER) || is.null(opts$PASS) || is.null(opts$URL)) stop("!!")

setwd(paste0(github,"/streamformatics/R"))


args = as.numeric(commandArgs(trailingOnly=TRUE))
if(length(args)==0){
  start = 1; end = 5000
}else{
  print(args)
  step = args[2]; 
  start = args[1]*step; end = start + step -1; 
}
dir.create("results")
rds_out =paste("results", paste("temp_depmap",start,end,"rds",sep="."),sep="/")



source("libs.R")
opts = setOptions("~/.sfx")
dbDir = opts$dbDir
keys = keyEnv$new(dbDir)
all = allEnv$new(dbDir,keys, "fspls/data.R")
uploadData=F
user = opts$USER

#depmap="/home/unimelb.edu.au/lcoin/Data/Depmap"

depmap="/data/gpfs/projects/punim1140/telsabeh/data_clean/data/depmap_data_primary"

#depmap="/home/unimelb.edu.au/lcoin/Data/Depmap"
dEnv = depmapEnv$new(depmap)

##UPLOAD DATA ONLY DO ONCE
#if(uploadData){
#  all_dist = allEnv$new(dbDir,keys,"dist.R")
#  dist1 = all_dist$get("Coin","Depmap","primary",reset=F)
#  dist1$sampleID()
#  dEnv$importData(dist1)
 # dEnv$importPheno(dist1)
#}
pows = seq(0.2,3,by=.4)
exp_y = c()  # c(0.5,1.5,2.0, 2.5,3.0)
exp_x = c()


flags0 = list(nrep=1,batch=0,bigmatrix=F,splitBy="disease",nphenos=10000, 
             start = 1,
             keep="pancreatic",var_thresh_x_quantile = 0.5, var_thresh_y_quantile = 0.5)

transform_y = getYTransform(pow = pows, offset=0.1, norm=1000, n_random=20, CHECK=T)
flags0$transform_y=toJSON(transform_y)
#flags0[['transform']] = '{"x" :"function(x) x", "exp":"function(x) exp(x)"}'
#flags0$transform=toJSON(getXTransform(c(-.5,.5,1.0,1.5)))
if(file.exists(rds_out)){
  rds1 = readRDS(rds_out)
}else{
rds1=readRDS("temp_depmap1.rds")
}
if(!is.null(rds1)) flags0=rds1$flags0

datasAll = all$get1("Coin","Depmap",flags=flags0, reload=T)
#sigDB = datasAll$getSigDB(db, user=user)
#data_flags = sigDB$get_data_flags()
#flags0 = fromJSON(data_flags$flags[[1]])
##reload
#datasAll = all$get1("Coin","Depmap",flags=flags0, reload=T)

#datasAll$updateTransforms(flags0$transform_y)
#pows = 1


phens_all = datasAll$pheno(sep=T,sep_group=F)
#phens_all = datasAll$pheno(sep=T)
#phens_all = dEnv$getPhens(datasAll, topphens)
#phens_all = datasAll$pheno(sep=F, memb=memb)
#phens_all = dEnv$getPhens(datasAll, topphens,sep=T)
options("fspls.types"=
          fromJSON('{"gaussian": ["correlation","var","mad","rms"],"binomial":"AUC","multinomial":"AUC","ordinal" : "AUC"}'))
#phens = phens[1:4]
flags = list(var_quantile = 0.00,# genes_incls =genes_incls,
              quantiles ="[0.0]" ,
              only_all=T,min=0,max=10,
              pthresh = 0.05, topn=20,useglmnet=T,stop_y="rand",
             nrep=0,batch=1, beam=1, train=names(datasAll$datas)) ## return can be model, vars or eval
 #vars = datasAll$convert(genes_incls,phens)
if(!is.null(rds1$flags)) flags = rds1$flags
vars_all1 = list(); all_models1 = list(); eval01= list(); ggps_l=list()


#sigDB$drop_all()
#area_p_l = vector(mode = "list", length = length(phens_all));
#plots = vector(mode = "list", length = length(phens_all));
#names(area_p) =unlist( lapply(phens_all, function(p1) names(p1[[1]])))
#names(plots) = names(area_p)
#sigDB$clear_all_user(user);  #dangerous  ..clears everything
#rds1 = readRDS("temp_depmap1.rds")
rds = rds1
if(!is.null(rds1) && file.exists(rds_out)){
 evals_all = rds$evals_all
 area_p_l = rds$area_p_l
 #start = rds$kk
}else{
#  start = 1
evals_all = list()
area_p_l = list()
}
db=paste("combined",start,sep=".")
#saveRDS(list(flags0=flags0, phens_all = phens_all, flags = flags),"flags.rds")
for(kk in start:min(end,length(phens_all))){
  print(paste(kk,"of", length(phens_all)))
 
  phens1 = phens_all[[kk]];phens = phens1
  print(phens1);
  title = phens1$gaussian[[1]]
  title = dEnv$convertFromCode(title)
if(!is.null(area_p_l[[title]])){
  print(paste("done",title)); next;
}else{
  print(paste("doing",title))
  
}

#f1=sigDB$flags(phens = phens1)
  #sigDB$clear_results(flags, phens1, transform_y) -- this will clear this result from sigDB
flags$nrep=1; flags$batch=0
flags$max=10
#sigDB$experiments(all_cols=T)
    vars_all = datasAll$select ( phens1, flags,verbose=F,db=db,user=user)
    vars_all21 = .extractFullVars(vars_all)
    
    if(length(vars_all21$variables)==0){
      print(paste("FOUND NOTHING", title,kk))
      next;
    }
    print(vars_all21$variables)
#    print(names(vars_al21$variables))
    flags$max = length(vars_all21$variables[[1]])
    flags$nrep=0; flags$batch=1
    vars_all = datasAll$select ( phens1, flags,verbose=F,db=db,user=user)
    
    #if(length(vars_all1$variables)==0) next;
    all_models =datasAll$makeAllModels(vars_all,flags=flags,verbose=F, db=NULL,user=user)
    eval1 = datasAll$evaluateAllModels(all_models,db=db,user=user)
    evals_all[[title]]=eval1
 #   ggps=.plotAll(evals_all, 36)
#    print(ggps[[1]])
    if(kk %% 10 ==0){
      print(paste("SAVING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!",kk))
    rds = list(evals_all = evals_all, area_p_l = area_p_l, kk=kk,    flags=flags, flags0=flags0)
    saveRDS(rds, rds_out)
    }
    if(FALSE){
      ggps=.plotEval1(eval1,legend=T, grid1=c("subpheno","pheno"), grid0="measure",
                       shape_color=c("cv_full"),
                      sep_by=c("data"),
                      showranges=T,linetype="cv_full",
                       scales="free",title =names(phens)[1], title1=title)
      ggps1[[title]]=ggps[[1]]
      print(ggps[[1]])
    }
    if(!is.null(eval1) && FALSE){
#      ev3=.takeAvg(eval1)
      
      
     ggps1=.plotEval2(eval1,legend=F, grid1="", grid0="measure",linetype="cv_full",
                     shape_color=c("pheno","data"),
                     #sep_by=c("cv_full"),
                     showranges=F, scales="free",title =title, title1="pheno" ) #, grid="pheno~cv_full",showranges = F)
#    ggps1[[1]]$`CV= FALSE FULL= TRUE`
      print(ggps1[[1]]$`CV=avg`)    
      ggps1[[1]]$`CV= FALSE FULL= TRUE`
    }
    if(TRUE){
    predictions_CV =datasAll$extractPredictions(all_models, CV = T, liab=F, data_nme = names(datasAll$datas)[[1]]);
    predictions =datasAll$extractPredictions(all_models, CV = F, liab=F, data_nme = names(datasAll$datas)[[1]]);
    
    area_p = .getAreaPlot1(predictions, families="gaussian") %>% tibble::add_column(CV="No CV", liab=F)
    
    if(length(predictions_CV)>0){
      area_CV = .getAreaPlot1(predictions_CV, families="gaussian") %>% tibble::add_column(CV="Leave-one-out CV", liab=F)
      area_p = rbind(area_CV, area_p)
    }
    if(FALSE){
      area_p$CV = factor(area_p$CV);area_p$lens = factor(area_p$lens)
      ggplot(area_p, aes_string(x="knots", y="value", color="lens"))+geom_point()+facet_grid("CV~lens")
    }
    area_p_l[[title]] = area_p
    #.compareCorrelation(area_p)
}
  
    #aa=roc(predictions[[2]]$y, predictions[[2]]$X0)
   
 
}

rds = list(evals_all = evals_all, area_p_l = area_p_l, kk=kk,    flags=flags, flags0=flags0)
saveRDS(rds, rds_out)

print("finished")
