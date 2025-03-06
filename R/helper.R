

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




.GET<-function(endpoint, org=NULL, project=NULL, db=NULL, query = list(),body = list(), outp="json", encode="json"){
  auth = authenticate(getOption("USER"), getOption("PASS"));baseurl = getOption("URL")
  url=paste(baseurl,endpoint,sep="/")
  if(!is.null(org)) url = paste(url, org,sep="/")
  if(!is.null(project)) url = paste(url,project,sep="/")
  if(!is.null(db)) url = paste(url,db,sep="/")
  print(url)
  #body = list()
  #if(length(flags)>0) body = list(flags = toJSON(flags))
  #print(body)
  a=GET(url, query = query,body=body,encode=encode,auth)
  a1=rawToChar(a$content)
  if(outp=="json") return(fromJSON(a1))
  if(outp=="tsv") return(readr::read_tsv(a1))
  return(a1)
}

.POST<-function(endpoint, org=NULL, project=NULL, db=NULL, query=list(),body = list(), files=NULL, outp="json",encode="multipart"){
  auth = authenticate(getOption("USER"), getOption("PASS"));baseurl = getOption("URL")
  url=paste(baseurl,endpoint,sep="/")
  if(!is.null(org)) url = paste(url, org,sep="/")
  if(!is.null(project)) url = paste(url,project,sep="/")
  if(!is.null(db)) url = paste(url,db,sep="/")
  print(url)
  #  body = list(flags=toJSON(flags))
  if(length(files)>0){
    if(length(files)>1){
      names(files) = rep("upload", length(files))
      uploads=lapply(files, function(f) upload_file(f))
      body = c(body, uploads)
    }else{
      body$upload =upload_file(files[[1]])
    }
  }
  if(length(query)>0){
    a=POST(url,body = body , query=query,auth,encode=encode) #,encode="multipart")
  }else{
    a = POST(url, body = body, query = query, auth, encode=encode)
  }
  a1=rawToChar(a$content)
  if(outp=="json") return(fromJSON(a1))
  if(outp=="tsv") return(readr::read_tsv(a1))
  return(a1)
}
