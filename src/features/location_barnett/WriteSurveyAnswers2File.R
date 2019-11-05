WriteSurveyAnswers2File <-
function(fildir,survey_id,SID){
  try1=try(setwd(fildir),silent=TRUE)
  if(class(try1) == "try-error"){
    warning(paste(fildir,"does not exist."))
    return(NULL)
  }
  filelist <- list.files(pattern = "\\.csv$")
  if(length(filelist)==0){return(NULL)}
  date_v=substr(filelist,1,10)
  ## keep only the last survey of each day
  IDkeep=length(date_v)
  curdate = date_v[length(date_v)]
  for(i in length(date_v):1){
    nexdate=date_v[i]
    if(curdate!=nexdate){
      IDkeep=c(IDkeep,i)
      curdate=nexdate
    }
  }
  IDkeep=rev(IDkeep)
  
  outmat=c()
  for(i in IDkeep){
    x=read.csv(filelist[i],fileEncoding="UTF-8")
    outmat=rbind(outmat,c(date_v[i],x$answer))
  }
  try1=try(colnames(outmat)<-c("Date",paste(as.character(x$question.text)," (",as.character(x$question.answer.options),")",sep="")),silent=TRUE)
  if(class(try1) == "try-error"){
    warning(paste("Survey non-constant over time for ID:",SID))
    return(NULL)
  }  
  write.table(outmat,paste("SurveyAnswers_",survey_id,"_",SID,".txt",sep=""),quote=F,col.names=T,row.names=F,sep="\t")
  return("success")
}
