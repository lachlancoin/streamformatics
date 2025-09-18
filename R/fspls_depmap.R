
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
library(readr); library(cowplot); library(ggrepel)

##SET OPTIONS
if(is.null(opts$USER) || is.null(opts$PASS) || is.null(opts$URL)) stop("!!")

setwd(paste0(github,"/streamformatics/R"))
source("libs.R")
opts = setOptions("~/.sfx")
dbDir = opts$dbDir
keys = keyEnv$new(dbDir)
all = allEnv$new(dbDir,keys, "fspls/data.R")
uploadData=F
user = opts$USER

#depmap="/home/unimelb.edu.au/lcoin/Data/Depmap"
depmap="/data/gpfs/projects/punim1140/telsabeh/data_clean/data/depmap_data_primary"

dEnv = depmapEnv$new(depmap)

##UPLOAD DATA ONLY DO ONCE
#if(uploadData){
#  all_dist = allEnv$new(dbDir,keys,"dist.R")
#  dist1 = all_dist$get("Coin","Depmap","primary",reset=F)
#  dist1$sampleID()
#  dEnv$importData(dist1)
 # dEnv$importPheno(dist1)
#}
  
flags0 = list(nrep=1,batch=0,bigmatrix=F,splitBy="disease",nphenos=10000, 
             start = 1,
             keep="lung",var_thresh_x_quantile = 0.5, var_thresh_y_quantile = 0.5)
#flags0[['transform']] = '{"x" :"function(x) x", "exp":"function(x) exp(x)"}'
t_x=getXTransform(c(-.2, .2,1.0,1.5))
flags0[['transform']] =toJSON(t_x)#  '{"x" :"function(x) x", "log":"function(x) log1p(x)"}'

datasAll = all$get1("Coin","Depmap",flags=flags0, reload=T)
datasAll$dims()

topphens = c("cyclovalone","LY2334737","astragaloside-a",
"10-hydroxycamptothecin","orotic-acid",
"maltobionic-acid","talmapimod","cyclopenthiazide",
"indirubin-3-monoxime","indoprofen","chlorhexidine",
"piromidic-acid","oleanolic-acid","epinephrine",
"homoquinolinic-acid","tipranavir","lidocaine","cyt387",
"GW-788388","midafotel"
)

#phens = datasAll$pheno(sep=F,sep_group=F)
phens_all = datasAll$pheno(sep=T)
phens_all = dEnv$getPhens(datasAll, topphens)
options("fspls.types"=
          fromJSON('{"gaussian": ["correlation","var","mad","rms"],"binomial":"AUC","multinomial":"AUC","ordinal" : "AUC"}'))
#phens = phens[1:4]
flags = list(var_quantile = 0.00,# genes_incls =genes_incls,
              quantiles ="[0.0]" ,
              only_all=T,min=0,max=10,
              pthresh = 1e-7, topn=20,useglmnet=T,
             nrep=10,batch=0, beam=1, train=names(datasAll$datas)) ## return can be model, vars or eval
#transform_y =  c("function(y) y","function(y) y")
#transform_y =  c("function(y) exp(-y)","function(y) -1*log(y)")
transform_y = getYTransform(pows = c(1),n_random=1)


#funcstrs = lapply(transform_y,function(x) eval(str2lang(x) ))
#transform_y =  c("function(y) exp(y)","function(y) log(y)")


#options("fspls.family"=NULL)
 #vars = datasAll$convert(genes_incls,phens)
vars_all1 = list(); all_models1 = list(); eval01= list(); ggps_l=list()
sigDB = datasAll$getSigDB("combined", user=user)

#area_p_l = vector(mode = "list", length = length(phens_all));
#plots = vector(mode = "list", length = length(phens_all));
#names(area_p) =unlist( lapply(phens_all, function(p1) names(p1[[1]])))
#names(plots) = names(area_p)
#sigDB$clear_all_user(user);  #dangerous  ..clears everything
for(kk in 1:length(phens_all)){
  print(paste(kk,"of", length(phens_all)))
  phens1 = phens_all[[kk]];phens = phens1
  print(phens1);
  title = names(phens1$gaussian)

#f1=sigDB$flags(phens = phens1)
  #sigDB$clear_results(flags, phens1, transform_y) -- this will clear this result from sigDB
    vars_all = datasAll$select ( phens1, flags,transform_y=transform_y,verbose=F,db="combined",user=user)
    #vars_all21 = .extractFullVars(vars_all2)
    #if(length(vars_all1$variables)==0) next;
    all_models =datasAll$makeAllModels(vars_all,verbose=F, db="combined",user=user)
    eval1 = datasAll$evaluateAllModels(all_models,db="combined",user=user)
    print(vars_all$variables)
    if(!is.null(eval1)){
     ggps1=.plotEval2(eval1,legend=T, grid1="", grid0="measure",linetype="pheno",shape_color=c("pheno","data"),sep_by=c("cv_full"), showranges=T, scales="free",title =title, title1="pheno" ) #, grid="pheno~cv_full",showranges = F)
#    ggps1[[1]]$`CV= FALSE FULL= TRUE`
      print(ggps1[[1]]$`CV=avg`)    
      #ggps1[[1]]$`CV= FALSE FULL= TRUE`
      #ggps1[[1]]$`CV= FALSE FULL= TRUE`
    }
    
   
  
    #aa=roc(predictions[[2]]$y, predictions[[2]]$X0)
   
 
}

## find goo results
evals=sigDB$evals(flags=flags,transform_y = transform_y)
evals1 = subset(evals, measure=="correlation")
evals2 = subset(evals1, cv_full=="CV=avg")
evals3 = subset(evals2,nsamps>5 & numvars>4 )
hist(evals3$mid)
evals4 = subset(evals3, mid>0.5)
topphens=dEnv$convertFromCode( unique(evals4$pheno))
match(evals4$pheno, unlist(phens_all))


print(topphens)
##
area_p_l = vector(mode = "list", length = length(phens_all));
phens_all_r = dEnv$getPhens(datasAll,topphens,sep=T)
## get sample correlation
for(kk in 1:length(phens_all)){
    print(paste(kk,"of", length(phens_all)))
    phens1 = phens_all[[kk]];phens = phens1
    print(phens1);
    
#    sigDB$evals(flags=flags, phens = phens1);
 #   mods=sigDB$loadModels(flags=flags, phens = phens1, transform_y = transform_y)
    
  #  vars_all = datasAll$select ( phens1, flags,transform_y=transform_y,verbose=F,db="combined",user=user)
    all_models =datasAll$makeAllModels(NULL, phens = phens1, flags = flags, transform_y = transform_y, verbose=F, db="combined",user=user)
    if(length(all_models$models)==0) next; ## just use phenotypes where we learned a model
    predictions_CV =datasAll$extractPredictions(all_models, CV = T, liab=F, data_nme = names(datasAll$datas)[[1]]);
    predictions =datasAll$extractPredictions(all_models, CV = F, liab=F, data_nme = names(datasAll$datas)[[1]]);
    
    area_p = .getAreaPlot1(predictions, families="gaussian") %>% tibble::add_column(CV="No CV", liab=F)
    if(length(predictions_CV)>0){
      area_CV = .getAreaPlot1(predictions_CV, families="gaussian") %>% tibble::add_column(CV="Leave-one-out CV", liab=F)
      area_p = rbind(area_CV, area_p)
    }
  
   # area_p$pheno = title
    area_p_l[[kk]] = area_p
    if(FALSE){
      ggp_pred0=.plotArea(area_p, rename=F,grid="CV~lens", CV=T,addText=T, scales="fixed", max_vars=3,r2=T)
      ggp_pred0<-ggp_pred0+xlab("Observed viability")+ylab("Predicted viability")+ theme(legend.position = legend_position,legend.title = element_text(size = txtsize))+ggtitle(title)
      print(ggp_pred0)
      filename=paste(title,"pdf",sep=".")
      ggsave(filename, ggp_pred0)
      plots[[kk]] = ggp_pred0
    }
}
area_p_l = area_p_l[unlist(lapply(area_p_l, length))>0]
names(area_p_l) = dEnv$convertFromCode(unlist(phens_all))
area_p_l_max =.merge1_new( lapply(area_p_l, .takeMax1, max_vars=100),addName="drug")
area_p_l_base =.merge1_new( lapply(area_p_l, .takeMax1, max_vars=0),addName="drug")
area_p_l_max$pheno = area_p_l_max$drug
area_p_l_base$pheno = area_p_l_base$drug
#rdsf=readRDS("area_p_l.rds"); area_p_l = rdsf$area_p_l; area_p_l_max = rdsf$area_p_l_max
#saveRDS(list(area_p_l=area_p_l, area_p_l_max=area_p_l_max, area_p_l_base=area_p_l_base), file="area_p_l.rds")
print(dim(area_p_l_max))
print(dim(area_p_l_base))
inds2 =which(area_p_l_max$lens>0)
area_p_l_max = area_p_l_max[inds2,,drop=F]
area_p_l_base = area_p_l_base[inds2,,drop=F]
#print(dim(area_p))
#area_p = .merge1_new(area_p_l, addName="pheno")
#area_p$pheno = factor(area_p$pheno)
#area_p$CV = factor(area_p$CV, levels = rev(unique(area_p$CV)))
#area_p0=area_p
legend_position="bottom"; txtsize=1;
options("ggrepel.max.overlaps"=1000)
max_vars = 10; max_samps=9; maxphens=5000
area_lens = list("base"=area_p_l_base,"fitted"=area_p_l_max)
nme_lens = names(area_lens); names(nme_lens)=nme_lens
nmel1 = nme_lens[[1]]
phens_list = names(area_p_l)
cv_l = unique(area_p_l2$CV); names(cv_l) = cv_l
cv_l1 = cv_l[[1]]
funcstr1 = eval(str2lang(transform_y[2]))
ggp_preds=lapply(cv_l, function(cv_l1){
plots2 = lapply(nme_lens, function(nmel1){
  area_p_l2 = area_lens[[nmel1]]
 
    area_p0 = subset(area_p_l2, CV==cv_l1 )
    area_p0$value= funcstr1(area_p0$value)
    title = paste(cv_l1, nmel1)
    print(title)
    print(cor(area_p0[,1:2],use="pairwise.complete.obs"))
    print(calcRMS(area_p0$knots, area_p0$value))
    print(calcVar(area_p0$knots, area_p0$value))
    dim(area_p0)
   # area_p0 = subset(area_p0, pheno %in% topphens)
  #  dim(area_p0)
  ggp_pred0=.plotArea(area_p0, arrange_by_sample=T,rename=F,grid="sample",, CV=T,addText=T,maxphens = maxphens,takeMax=F,showText=F,reorder=T,
                      code_len=3,shapes=F, scales="fixed", max_vars=max_vars, max_samps = max_samps,r2=T,p=c(0.01,0.99))
  ggp_pred0+xlab("Observed viability")+ylab("Predicted viability")+ theme(strip.text = element_text(size = 10),legend.position = legend_position,legend.title = element_text(size = txtsize))+ggtitle(title)
})
ggp1 = plot_grid(plots2[[1]], plots2[[2]])
filename=paste("combined_sample_plot",max_vars,max_samps, maxphens,cv_l1, "pdf",sep=".")
ggsave(filename, ggp1,  width=100, height=200, units="cm",limitsize=F)
})
#print(ggp_preds$`No CV`)

    




