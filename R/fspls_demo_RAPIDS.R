
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

opts = setOptions("~/.sfx")
rlib = opts$rlib
dbDir = opts$dbDir
github=opts$github
.libPaths(rlib)

options(bigmemory.allow.dimnames=TRUE)
#install.packages("readr",lib=rlib)
library(curl);library(httr); library(jsonlite);library(ggplot2)
library(readr);

##SET OPTIONS
if(is.null(opts$USER) || is.null(opts$PASS) || is.null(opts$URL)) stop("!!")

setwd(paste0(github,"/streamformatics/R"))
source("libs.R")
opts = setOptions("~/.sfx")
dbDir = opts$dbDir
keys = keyEnv$new(dbDir)
all = allEnv$new(dbDir,keys, "fspls/data.R")

### MAKE NEW DATABAES
all_dist = allEnv$new(dbDir,keys,"dist.R")
keys$makeDB("lcoin","Coin","LCAH","Discovery")



dist1 = all_dist$get("Coin","LCAH","Discovery")
dist2 = all_dist$get("Coin","LCAH","Validation")
flags0 = list(nrep=1, bigmatrix=T, phen_nme="all",all_v_all=T,merge=T, 
              duplicate_ordinal=c("gaussian"),
              var_thresh_x_quantile=0.5) #"gn"
old_sigs = read_json("~/Data/sAPI/disease_severity_signature_weights.json")
#old_sigs = read_json("~/Data/sAPI/disease_class_signature_weights.json")
genes_incls= unlist(lapply(old_sigs$weights$Sig.1$attributes$dimnames$variables,function(str)strsplit(str,"_")[[1]][1]))[-1]
#flags[['transform']] =toJSON(old_sigs$attributes$normalise)
flags0$genes_incls = toJSON(genes_incls)
steps=c(seq(-1,2,by=0.4)) #,seq(1,10), seq(1,2,by=.5))
steps = c(1)

flags0[['transform']] =toJSON(.getTransformFuncs(steps))
flags0$transform=NULL
#flags0[['transform']] = '{"x" :"function(x) x", "log":"function(x) log1p(x)"}'

#"ans":"function(x) {x1=(x+3/8); 2*sign(x1)*abs(x1)^(0.5)}",,"bin":"function(x){ x[x>1e6]=1e6; sqrt(1e6)*asin(sqrt(x/1e6))}"}'
#flags0[['genes_incls']]=toJSON(genes_incls) #toJSON(list(first=genes_incls,second="all"))
#flags0$genes_incls=NULL
datasAll=all$get1("Coin","LCAH",flags=flags0,reload=T)
print(datasAll$pheno())
print(datasAll$dims())
#datas2=all$get1("Coin","LCAH","Validation",flags=flags)
dir="/home/unimelb.edu.au/lcoin/Data/sAPI/Coin/LCAH/LCAH_data"
dir1="/home/unimelb.edu.au/lcoin/Data/sAPI/Coin/LCAH"
nme="22022024_RAPIDS_Pheonix_scores.csv"
nme1 = "Discovery_meta.csv"
nme2 = "Validation_meta.csv"
filenme1 = paste(dir1,nme1,sep="/")
filenme2 = paste(dir1,nme2,sep="/")

phenos = dist1$getPheno("all")
phenos2 = dist2$getPheno("all")
phenos$upload_pheno(filenme1,nme1,opts$USER,flags, samples=dist1$sampleID())
phenos2$upload_pheno(filenme2,nme2,opts$USER,flags, samples=dist2$sampleID())

phenos$getPhenos()
diffs = list(pheonix4=c("pheonix4_00","pheonix4_24"),pheonix8=c("pheonix8_00","pheonix8_24"))

#diffs = list(gn_od=c("gn_od_00","gn_od_24"))
phenos$importDiff(user,diffs)
phenos2$importDiff(user,diffs)

phenos$upload_pheno(filenme1,nme1,opts$USER,flags, samples=samples)
export = phenos$exportPheno()
#dist1$upload_pheno()
#####


#flags = list(nrep=10,batch=0)
#datasAll=all$get1("Coin","golub","golub_data",flags=flags)
#datasAll=all$get1("Coin","LCAH",flags=flags, reload=T)

##RUN
#datasAll=all$get1("Coin","LCAH",flags=flags,reload=T)

datasAll$dims()
user=opts$USER


phens = datasAll$pheno(sep=F,sep_group=F)
#phens = phens[[1]]
#phens2 = phens
to_rem = unique(c(grep("Probable",(phens2$binomial.multiway$binomial.multiway)),
grep("Combined",(phens2$binomial.multiway$binomial.multiway)),
grep("Unknown",(phens2$binomial.multiway$binomial.multiway))))
phens2$binomial.multiway$binomial.multiway = phens2$binomial.multiway$binomial.multiway[-to_rem]

phens  = phens$all[-7]
phens$gaussian = phens$gaussian[-7]

#phens = phens[3:6]
#phens$all = phens$all[1]
#phens = phens[grep("binomial.1", names(phens))]
#grep(bestpheno,unlist(phens))
#phens$all = phens$all[1]
#phens = phens[grep("multiway", names(phens))]
options("fspls.types"=
          fromJSON('{"gaussian": ["correlation","var","mad"],"binomial":"AUC","multinomial":"AUC","ordinal" : "AUC"}'))
#phens = phens[1:4]
flags1 = list(quantiles = toJSON(0.5),# genes_incls =genes_incls,
              project=T, useoffset=T, useglmnet=T,
              stop_y="rand",
              pthresh = 0.05,max=50, topn=1,nrep=10,batch=0, beam=1, train=names(datasAll$datas)[1]) ## return can be model, vars or eval
#flags1[['transform_y']]  = '{"y":"function(y) y","expy" : "function(y) exp(y)"}'
#flags1[['transform_y_inverse']]  = '{"y":"function(y) y","expy" : "function(y) log(y)"}'

#flags1$test=flags1$train
#datasAll$update(flags1)
#aa = fromJSON('{"rna_star_log.AATBC;rna_star_log.MAFG;rna_star_log.VAV1;rna_star_x.MS4A7;rna_star_log.MPP7;rna_star_log.ATP6V0A1":{"rna_star_log.AATBC":["rna_star_log","AATBC"],"rna_star_log.MAFG":["rna_star_log","MAFG"],"rna_star_log.VAV1":["rna_star_log","VAV1"],"rna_star_x.MS4A7":["rna_star_x","MS4A7"],"rna_star_log.MPP7":["rna_star_log","MPP7"],"rna_star_log.ATP6V0A1":["rna_star_log","ATP6V0A1"]}}')
#genes_incls = aa
options("fspls.verbose1"=T);options("fspls.check"=T)
options("fspls.family"=NULL)
 #vars = datasAll$convert(genes_incls,phens)
transform_y=getYTransform(n_random=5)

#transform_y = list(x=c("function(y) y", "function(y) y"))
#phens1=phens$all
sigDB = datasAll$getSigDB("combined")
#sigDB$drop_all()
options("fspls.check"=T)
#phens$gaussian = phens$gaussian[1]
  vars_all = datasAll$select ( phens, flags1,transform_y, verbose=T,db="combined")
  warnings()
  #options("fspls.family"=NULL)
  vars_all_full=.extractFullVars(vars_all); print(vars_all_full$variables)
#  options("glmnet"=T)
  all_models =datasAll$makeAllModels(vars_all)
  #all_models1 =datasAll$makeAllModels(vars,phens,flags1)
  
  flags1[['liab']]=T
 # eval0 = datasAll$evaluateAllModels(all_models,verbose=F)
#  eval1 = datasAll$evaluateAllModels(all_models,verbose=F)
  eval3 = datasAll$evaluateAllModels(all_models,verbose=F)
  print(eval3$experiment_id)
  #comb=rbind(eval0,eval1,eval2,eval3)
  #comb$experiment_id = factor(comb$experiment_id)
  
 
  
  #.calcAverageAccuracy(comb)
  ggps5=.plotEval1(eval3,legend=T, grid0="pheno", grid1="measure", point=T,line=F,shape_color=c("experiment_id", "data", "subpheno"),
                   sep_by=c("cv_full"), showranges=F, scales="free") #, grid="pheno~cv_full",showranges = F)
ggps5$`CV=avg`
ggps5$`CV= FALSE FULL= TRUE`

  plot_grid(ggps5$`CV=avg 0`, ggps5$`CV=avg 1`)
  
  
   # subset(eval0, isfull &cv)
  
  #toSave = list(vars = vars, all_models = all_models, eval0 = eval0)  
  #saveRDS(toSave,"all1.rds")
  eval_cv_full = subset(eval0,isfull & cv & measure=="correlation" ) 
  eval_cv_full = subset(eval0,isfull & cv & measure=="correlation" ) 
  
  hist(eval_cv_full$mid,br=100)
  eval_cv = (subset(eval0,cv==F & isfull & numvars==2 & measure=="correlation"))
 # which.max(eval_cv_full$mid)
 eval0=subset(eval0, measure=="correlation")
   eval1 = .calcEval1(eval0,rename=F)
  #eval1_avg = subset(eval1, cv=="CV= avg" & measure=="correlation")
   transf = unique(eval0$transform_y); names(transf) = transf
   ggps = lapply(transf, function(tf){
     .plotEval1(subset(eval1, transform_y==tf),legend=T, grid0="pheno", shape_color=c("trainedOn","data","measure"),showranges=F,
                  linetype=c("fullmodel","transform_y"),txtsize=0.1,scales="fixed",point=F,
                  sep_by=c("cv_full")) #, grid="pheno~cv_full",showranges = F)
   })
  ggps=.plotEval1(eval1,legend=T, grid="pheno~cv_full", shape_color=c("measure"),sep_by="", showranges=F) #, grid="pheno~cv_full",showranges = F)
    ggps=.plotEval1(eval1,legend=T, grid="subpheno~data", shape_color=c("measure"),sep_by=c("cv_full"), showranges=F) #, grid="pheno~cv_full",showranges = F)

    ggps5=.plotEval1(eval1,legend=T, grid0="pheno", grid1="subpheno", point=F,shape_color=c("data"),sep_by=c("cv_full"), showranges=F, scales="free") #, grid="pheno~cv_full",showranges = F)
    ggps5=.plotEval1(eval1,legend=T, grid1="", grid0="pheno", point=F,shape_color=c("data"),sep_by=c("cv_full"), showranges=F, scales="fixed") #, grid="pheno~cv_full",showranges = F)
    
  plot_grid(ggps1[[1]], ggps1[[2]], ggps1[[3]], ggps1[[4]])  
  plot_grid(ggps$`CV= FALSE FULL= TRUE`, ggps1$`CV= FALSE FULL= TRUE`,align="v",axis="v")
    
  ggps
  out = paste0("plot_species_",max(eval1$numvars),".pdf");
  pdf(out,width=30,height=30)
  for(ggp in ggps) print(ggp)
  dev.off()
  
  #eval1 = eval1[grep("avg", eval1$cv, inv=T),]
  eval1 = subset(eval1, cv_full=="CV= TRUE FULL= TRUE")
  ggp=.plotEval1(eval1, grid="measure~cv_full", rename=T,shape_color="pheno", linetype="fullmodel",logy=F, legend=F,showranges=F)
  eval1
#})
eval0 = .merge1_new(eval_all)
eval2 = subset(eval0,cv=="CV= avg" & measure=="var")
#eval2 = subset(eval_,cv=="CV= avg" & measure=="var")
topmid =tail(sort(eval2$mid,decr=T),1)
best_pheno=unique(eval2$pheno[match(topmid,eval2$mid)])

eval3 = subset(eval0,pheno %in% best_pheno)
eval3$pheno=factor(eval3$pheno, levels = best_pheno)
#eval2 = subset(eval,pheno==eval$p)


ggps=.plotEval1(eval3, grid="pheno~cv", rename=T,shape_color="measure", linetype="fullmodel",logy=T)
ggps
#ggps=.plotEval1(eval,sep="pheno", grid="data~cv", rename=T,shape_color="subpheno", linetype="fullmodel")

predictions0 =datasAll$extractPredictions(all_models,phens2, flags, CV = F, liab=F, data_nme = names(datasAll$datas)[[1]]);
predictions1 =datasAll$extractPredictions(all_models,phens, flags, CV = F, liab=F,data_nme = names(datasAll$datas)[[2]]);
#aa=roc(predictions[[2]]$y, predictions[[2]]$X0)
ggp_pred0=.plotArea1(predictions0, rename=T,max_vars=44)
ggp_pred0

ggp_pred1=.plotArea1(predictions1, rename=T,max_vars=44)
ggp_pred1
#.plotArea(predictions$y$all, rename=T)
ggps1 = list(ggp_pred0, ggp_pred1)
out = paste0("plot_predictions_",max(eval1$numvars),".pdf");
pdf(out,width=30,height=30)
for(ggp in ggps1) print(ggp)
dev.off()


datasAll$angles(vars,phens,flags)

#pseudorandom
rotl <- function(x, n, bit_length = 32) {
  if (n < 0 || n >= bit_length) {
    stop("Shift amount 'n' must be between 0 and bit_length - 1.")
  }
  
  # Calculate the rotated value
  rotated_value <- bitOr(bitShiftL(x, n), bitShiftR(x, bit_length - n))
  return(rotated_value)
}
rotr <- function(x, n, bit_length = 32) {
  if (n < 0 || n >= bit_length) {
    stop("Shift amount 'n' must be between 0 and bit_length - 1.")
  }
  
  # Calculate the rotated value
  rotated_value <- bitOr(bitShiftR(x, n), bitShiftL(x, bit_length - n))
  return(rotated_value)
}

num <- 1:10 # Binary: 0000...1010
rotated_num_left <- rotl(num, 32)
print(rotated_num_left) # Expected: 40 (Binary: 0000...101000)
rotr(rotated_num_left,2)
