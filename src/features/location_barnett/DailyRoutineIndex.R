DailyRoutineIndex <-
function(indday,mobmat,subsetinds_v,subsetdayofweek_v,subsetstarttime_v,tz,CENTERRAD){
  submat=matrix(mobmat[subsetinds_v[[indday]],],ncol=7)
  daydist_v = rep(NA,length(subsetinds_v))
  for(i in 1:length(subsetinds_v)){
    if(i == indday){next}
    daydist_v[i]=DayDist(indday,i,mobmat,subsetinds_v,subsetstarttime_v,tz,CENTERRAD)
  }
  if(length(which(!is.na(daydist_v)))==0){
    circscore=NA
  }else{
    circscore = mean(daydist_v,na.rm=T)    
  }
  if(subsetdayofweek_v[indday] == "Saturday" || subsetdayofweek_v[indday] == "Sunday"){
    IDcompare = c(which(subsetdayofweek_v=="Saturday"),which(subsetdayofweek_v=="Sunday"))
  }else{
    IDcompare = c(which(subsetdayofweek_v=="Monday"),which(subsetdayofweek_v=="Tuesday"),which(subsetdayofweek_v=="Wednesday"),which(subsetdayofweek_v=="Thursday"),which(subsetdayofweek_v=="Friday"))
  }
  if(length(IDcompare)==0 || length(which(!is.na(daydist_v[IDcompare])))==0){
    wkscore=NA
  }else{
    wkscore = mean(daydist_v[IDcompare],na.rm=T)    
  }
  return(list('cscore'=circscore,'wscore'=wkscore))
}
