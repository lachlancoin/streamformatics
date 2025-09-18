
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
library(readr); library(cowplot); library(tidyr)


setwd(paste0(github,"/streamformatics/R"))
source("libs.R")
opts = setOptions("~/.sfx")
dbDir = opts$dbDir
keys = keyEnv$new(dbDir)
org="Coin"; project="NITRIC"
keys$make_admin(opts$USER,org)

keys$makeDB(opts$USER,org,project,"Discovery")
all = allEnv$new(dbDir,keys, "fspls/data.R")


##LOAD DATA
all_dist = allEnv$new(dbDir,keys,"dist.R")
dirn="/home/lcoin/punim1140/Jingni/Batch_correction/DATA"
nmef = "NITRIC_Metadata_09112023.csv";
f= paste(dirn, nmef,sep="/")
#flags$id="record_id"
pheno =data.frame(data.table::fread(f, header=T, sep=","))
cohorts = unique(pheno$cohort); names(cohorts) = cohorts
sampletime = unique(pheno$Sampling_time); names(sampletime)=sampletime

filen = "NITRIC_corrected_counts.csv.gz"
mat = data.frame(data.table::fread(paste(dirn,filen,sep="/"),sep=","))
rownames(mat) = mat[,1]
colnames(mat) = unlist(lapply(colnames(mat), function(str) strsplit(str,"\\..")[[1]][1]))
mat1 = t(mat[,-1])
mat=mat1
flags = list()
dist_all = lapply(cohorts, function(cohort1){
  flags$id="record_id"
  keys$makeDB(opts$USER,org, project, cohort1)
  dist1 = all_dist$get(org,project,db=cohort1,reload=T, reset=T)
   lapply(sampletime, function(st1){
     print(st1)
    pheno_sub =  subset(pheno, cohort==cohort1 & Sampling_time==st1 )
   # phen_nme=paste(cohort1, st1, sep="::")
    st2 = tolower(st1)
    phenoE=dist1$getPheno(st2, refresh=T, clear=T); phenoE$clear_all();
    phenoE$importPheno(pheno_sub, nmef, opts$USER, flags)
    #aa = phenoE$readMatrix(p1)
    samps= phenoE$readMatrix("Sequencing_Sample_ID",replace_levels=T)
    samples = samps[,1]
    names(samples) = rownames(samps)
    res1 = dist1$importMat(mat, paste(filen,st1,sep="."), opts$USER, flags, type=st2, samples = samples)
    samp2 = dist1$sampleID()
    print(paste("here", length(samples),length(dist1$sampleID())))
    print(res1)
  #  p1=phenoE$getPhenos()
   # p2=dist1$readMatrix(st2)
  #  p2$all[1:5,1:5]
  #  grep("record_id", names(p1));
  #  grep("record_id", names(pheno_sub))
  #  phenoE$upload_pheno(pheno_sub,nmef,opts$USER,flags)
   })
   print(dist1$types())
   dist1
})

dist_all$validation$types()
s2=dist_all$discovery$phenos$pre_bypass$samples()[,2]
s1 = dist_all$discovery$sampleID()
##getPhens
phens_all=lapply(dist_all, function(d_all){
  lapply(sampletime, function(st1){
    st2 = tolower(st1)
    d_all$getPheno(st2, refresh=T)
  })
})
tbls = lapply(phens_all, function(p_all){
  lapply(p_all, function(p_all1){
    p_all1$getTable(NULL)
  })
})

dist_all$discovery$types()
##IMPORT DIFFS
deltas = lapply(dist_all, function(d_all){
  types = d_all$types()
  types1 = types$type
  names(types1) = types1
  tbls1 = lapply(types1, function(type){
      d1 =  d_all$readMatrix(type)
      d1$all
  })
  tbl1 = tbls1[[1]]
  dn1 = dimnames(tbl1)
  diffs=lapply(2:length(tbls1), function(k){
    tbl2 = tbls1[[k]]
    dn2 = dimnames(tbl2)
    mi2 = match(dn2[[2]],dn1[[2]])
    mi1 = match(dn2[[1]],dn1[[1]])    
    cbind(dn1[[1]][mi1[!is.na(mi1)]], dn2[[1]][!is.na(mi1)])
    diff = tbl1[mi1[!is.na(mi1)], mi2[!is.na(mi2)]] - tbl2[!is.na(mi1), !is.na(mi2)]
    res1 = d_all$importMat(diff, paste(filen,"diff",k,sep="."), opts$USER, flags, type=paste("diff",k,sep="."))
    res1
  })
  names(diffs) = names(tbls1)[-1]
  d_all$types()
  diffs
})

#flags1[['transform_y_inverse']]  = '{"expy" : "function(y) log(y)"}'
flags0 = list(nrep=1, bigmatrix=F, phen_nme=c("post_bypass"),all_v_all=F, max_ordinal=30,merge=T,
              duplicate_ordinal=c("binomial","gaussian"),
              var_thresh_x_quantile=0.5)
flags0[['transform']] = '{"x" :"function(x) x", "log1p":"function(x) log1p(x)"}'
#flags0[['transform']] =.getTransformFuncs(c(seq(-1,-0.2,by=0.2),1/seq(1,10), seq(1,2,by=.5)))#  '{"x" :"function(x) x", "log":"function(x) log1p(x)"}'


org="Coin"; project="NITRIC"
#old_sigs = read_json("/home/lcoin/disease_severity_signature_weights.json")
#genes_incls1= toJSON(unlist(lapply(old_sigs$weights$Sig.1$attributes$dimnames$variables,function(str)strsplit(str,"_")[[1]][1]))[-1])
#flags0$genes_incls = genes_incls
datasAll=all$get1("Coin", "NITRIC", flags=flags0, reload=T)
print(dimnames(datasAll$datas[[1]]$y$binomial)[[2]])

datasAll$dims()
quantile(datasAll$datas$discovery.all$vars$pre_bypass_x)


p1 = "perfusion_xclamp_total"; fam="gaussian"
p1="picu_pelod2_48"; fam=paste("ordinal",p1,sep=".")
p1="picu_pelod2_0"; fam=paste("ordinal",p1,sep=".")
p1="death_28day"; fam="binomial"
p1="vfd"; fam=paste("ordinal",p1,sep=".") 
p1="dem_gender"; fam="binomial"
p1="los_picu"; fam="ordinal"  
p1="mment_vent_dur"; fam="ordinal"  
p1="comp_outcome" ; fam="biomial"

p1 = c("death_28day","comp_outcome","vfd","los_picu","picu_lcos_any")

phens = datasAll$pheno(sep=F)
i1=unlist(lapply(p1, function(p2)grep(p2,phens$all)))
fam = unique(names(phens$all)[i1])
fam = fam[grep("ordinal",fam,inv=F)]  ## remove ordinal
phens$all =phens$all[names(phens$all) %in% fam]
for(j in 1:length(fam)){
  i2=unlist(lapply(p1, function(p2)grep(p2,phens$all[[fam[[j]]]])))
  phens$all[[fam[[j]]]] = phens$all[[fam[[j]]]][i2]
  
}
phens

options("fspls.family"=NULL)
options("fspls.types"=
          fromJSON('{"gaussian": ["correlation","var","mad","rms"],"binomial":"AUC","multinomial":"AUC","ordinal" : "AUC"}'))
flags1 = list(
              quantiles ="[0.0]" , ## remove bottom 10% variable
              only_all=T,pheno_balance=T,reweight=T,
              pthresh = 1e-5,max=1, topn=20,nrep=1,batch=0, beam=1, train=names(datasAll$datas)[1]) ## return can be model, vars or eval
flags1$data_types=toJSON(list("all"=names(datasAll$datas[[1]]$data)))
flags1$data_types = toJSON(list("post"=grep("post",names(datasAll$datas[[1]]$data),v=T)))
flags1$data_types = toJSON(list("pre"=grep("pre",names(datasAll$datas[[1]]$data),v=T),"all"=names(datasAll$datas[[1]]$data)))

phens1 = phens[[1]]
sigDB = datasAll$getSigDB("combined",reload=T)
flags2 = sigDB$flags(phens = phens1, flags = flags1, transform_y = transform_y)
if(FALSE && !is.null(flags2)){
    flags = flags2
    diff=.diffFlags(flags2, flags1)
    print(diff)
    flags1=flags2
    if(FALSE) sigDB$updateFlags(flags, flags1) ## updates the flags
}

sigDB = datasAll$getSigDB("combined")
#datasAll$update(phens1,flags1,  transform_y=c("function(y) y","function(y) y"))
#sigDB$clear_results(flags1, phens1, transform_y)
vars_all = datasAll$select ( phens1, flags1,verbose=T, db="combined")
all_models =datasAll$makeAllModels(vars_all, verbose=T, db="combined")
eval1 = datasAll$evaluateAllModels(all_models, db="combined")

eval1
vars_all1 = .extractFullVars(vars_all)
plots_sep=datasAll$plotData(vars_all1, phens1 = phens1, all_types=F, transform_x = flags0$transform, violin=F, assoc=F)
plots_sep
grid1="subpheno"
ggps1=.plotEval2(eval1,legend=T, grid1="", grid0="measure",linetype="pheno",shape_color=c("pheno","data"),sep_by=c("cv_full"), showranges=T, scales="free",title =names(phens)[1], title1="pheno" ) #, grid="pheno~cv_full",showranges = F)

ggps1[[1]]$`CV=avg`
ggps1[[1]]$`CV= FALSE FULL= TRUE`


all_models1 = .extractFullModels(all_models)

predictions0 =datasAll$extractPredictions(all_models, CV = F, liab=F, data_nme = names(datasAll$datas)[[1]]);
area_p = .getAreaPlot1(predictions0,families="binomial")
#aa=roc(predictions[[2]]$y, predictions[[2]]$X0)
ggp_pred0=.plotAreaSep(area_p, rename=F,sep="pheno")

ggp_pred0
all_sig = allEnv$new(dbDir,keys,"sigEnv.R")




