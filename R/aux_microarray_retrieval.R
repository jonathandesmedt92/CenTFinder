#######################################################################
###   CenTFinder                                                    ###
###   Auxiliary function definitions: Microarray retrieval          ###
###   Author: Jonathan De Smedt                                     ###
#######################################################################

# Include
#' @include classes.R

retrieveURLdata<- function(URL,FUN) {
  result <- tryCatch(FUN(url = URL), error=function(e){cat('In error handler\n'); print(e); e})
  return(result)
}

downloadURLdata<- function(URL, destfile) {
  tryCatch(download.file(url = URL, destfile = destfile), error=function(e){cat('In error handler\n'); print(e); e})
}

addFileNames<-function(sdrf){
  sdrf$FileNames<-apply(sdrf, MARGIN = 1, FUN = function(x){
    x<-as.character(x)
    x<-x[grepl(pattern = "\\.CEL",ignore.case = T,x=x)]
    if(length(x)==1){
      return(x)
    } else {
      return(NA)
    }
  })
  return(sdrf)
}

getSDRF<-function(AEcode){
  # Download
  sdrf<-tryCatch(data.table::fread(paste("ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/",substr(AEcode,3,6),"/",AEcode,"/",AEcode,".sdrf.txt",quote="",sep = "")),
                 error = function(e){NULL})
  if(is.null(sdrf)){
    return(NULL)
  }
  if(!nrow(sdrf)>0){
    return(NULL)
  }
  # Reformat
  sdrf<-data.frame(sdrf, stringsAsFactors = F)
  sdrf[]<-lapply(sdrf, as.character)
  names(sdrf)<-gsub(pattern = "\\.",replacement="",x=names(sdrf))
  sdrf<-addFileNames(sdrf)
  names(sdrf)<-tolower(names(sdrf))
  #sdrf<-sdrf[,names(sdrf)%in%c("assayname","characteristicsorganism","materialtype","protocolref","protocolref1","protocolref2","protocolref3","sourcename","technologytype","characteristicsorganismpart","characteristicscelltype","characteristicscellline","label","performer","dateofscanning","description","filenames")]
  sdrf$AEcode<-AEcode
  # Return
  return(sdrf)
}


