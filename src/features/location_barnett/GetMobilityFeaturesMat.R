GetMobilityFeaturesMat <-
function(mobmat,obj,mobmatmiss,tz,CENTERRAD,ITRVL){
  #### Get significant locations
  slout=SigLocs(mobmat,obj,CENTERRAD,tz=tz)
  IDhome=which(slout[,4]==1)
  if(length(IDhome)==0){IDhome=1}
  homex=slout[IDhome,1];homey=slout[IDhome,2]
  #### Partion mobmat data into daily subsets: create subsetinds_v and daystr_v.
  curdate=strsplit(as.character(as.POSIXct(mobmat[1,4],tz=tz,origin="1970-01-01"))," ")[[1]][1]
  curtime=strsplit(strsplit(as.character(as.POSIXct(mobmat[1,4],tz=tz,origin="1970-01-01"))," ")[[1]][2],":")
  subsetinds_v = list()
  subsetdayofweek_v = c()
  subsetstarttime_v = c()
  daystr_v = c(curdate)
  dayind=1
  subsetinds = c(1)
  if(nrow(mobmat)<2){
    outmat=NULL
    return(list(outmat,slout))
  }
  for(i in 2:nrow(mobmat)){
    nexdate=strsplit(as.character(as.POSIXct(mobmat[i,4],tz=tz,origin="1970-01-01"))," ")[[1]][1]
    if(nexdate == curdate && i<nrow(mobmat)){
      subsetinds=c(subsetinds,i)
    }else{
      subsetinds_v[[dayind]]=subsetinds
      subsetdayofweek_v = c(subsetdayofweek_v,weekdays(as.Date(curdate)))
      subsetstarttime_v = c(subsetstarttime_v,as.numeric(as.POSIXct(paste(curdate," 00:00:00",sep=""),tz=tz,origin="1970-01-01")))
      dayind=dayind+1          
      if(mobmat[i-1,1]==2 && (mobmat[i-1,7]-mobmat[i-1,4])/(60*60*24)>1){
        ii=1
        middate = strsplit(as.character(as.POSIXct(mobmat[i-1,4]+(60*60*24)*ii,tz=tz,origin="1970-01-01"))," ")[[1]][1]
        while(middate!=nexdate){
          subsetdayofweek_v = c(subsetdayofweek_v,weekdays(as.Date(middate)))
          subsetstarttime_v = c(subsetstarttime_v,as.numeric(as.POSIXct(paste(middate," 00:00:00",sep=""),tz=tz,origin="1970-01-01")))
          subsetinds_v[[dayind]]=i-1
          dayind=dayind+1
          if(middate!=daystr_v[length(daystr_v)]){
            daystr_v=c(daystr_v,middate)        
          }          
          ii=ii+1
          middate = strsplit(as.character(as.POSIXct(mobmat[i-1,4]+(60*60*24)*ii,tz=tz,origin="1970-01-01"))," ")[[1]][1]
        }
      }
      curdate=nexdate
      subsetinds=c(i-1,i)
      if(curdate!=daystr_v[length(daystr_v)]&& !(length(daystr_v)==length(subsetinds_v) && nrow(mobmat)==i)){
        daystr_v=c(daystr_v,curdate)        
      }
    }
  }
  #### Partion mobmatmiss data into daily subsets: create subsetinds_v and daystr_v.
  curdate=strsplit(as.character(as.POSIXct(mobmatmiss[1,4],tz=tz,origin="1970-01-01"))," ")[[1]][1]
  curtime=strsplit(strsplit(as.character(as.POSIXct(mobmatmiss[1,4],tz=tz,origin="1970-01-01"))," ")[[1]][2],":")
  subsetindsmiss_v = list()
  daystrmiss_v=c(curdate)
  dayind=1
  subsetinds = c(1)
  for(i in 2:nrow(mobmatmiss)){
    nexdate=strsplit(as.character(as.POSIXct(mobmatmiss[i,4],tz=tz,origin="1970-01-01"))," ")[[1]][1]
    if(nexdate == curdate && i<nrow(mobmatmiss)){
      subsetinds=c(subsetinds,i)
    }else{
      curdate=nexdate
      subsetindsmiss_v[[dayind]]=subsetinds
      if(curdate!=daystrmiss_v[length(daystrmiss_v)] && !(length(daystrmiss_v)==length(subsetindsmiss_v) && nrow(mobmatmiss)==i)){
        daystrmiss_v=c(daystrmiss_v,curdate)        
      }
      subsetinds=c(i-1,i)
      dayind=dayind+1
      if(mobmatmiss[i-1,1]==2 && (mobmatmiss[i-1,7]-mobmatmiss[i-1,4])/(60*60*24)>1){
        ii=1
        middate = strsplit(as.character(as.POSIXct(mobmatmiss[i-1,4]+(60*60*24)*ii,tz=tz,origin="1970-01-01"))," ")[[1]][1]
        while(middate!=nexdate){
          subsetindsmiss_v[[dayind]]=i-1
          dayind=dayind+1
          if(middate!=daystrmiss_v[length(daystrmiss_v)]){
            daystrmiss_v=c(daystrmiss_v,middate)        
          }          
          ii=ii+1
          middate = strsplit(as.character(as.POSIXct(mobmatmiss[i-1,4]+(60*60*24)*ii,tz=tz,origin="1970-01-01"))," ")[[1]][1]
        }
      }
    }
  }
  ##### intersect mobmat and mobmatmiss to ignore missing data
  IDkeep=c()
  IDkeepmiss=c()
  for(i in 1:length(daystr_v)){
    for(j in 1:length(daystrmiss_v)){
      if(daystr_v[i]==daystrmiss_v[j]){
        IDkeep=c(IDkeep,i)
        IDkeepmiss=c(IDkeepmiss,j)
        break
      }
    }
  }
  if(length(IDkeep)==0){
    outmat=NULL
    return(list(outmat,slout))
  }
  daystr_v=daystr_v[IDkeep]
  subsetinds_v=subsetinds_v[IDkeep]
  subsetdayofweek_v=subsetdayofweek_v[IDkeep]
  subsetstarttime_v=subsetstarttime_v[IDkeep]
  daystrmiss_v=daystrmiss_v[IDkeepmiss]
  subsetindsmiss_v=subsetindsmiss_v[IDkeepmiss]
  ##### Compute mobility features for each day
  Nfeatures=15
  outmat=matrix(,nrow=length(daystr_v),ncol=Nfeatures)
  colnames(outmat)=c("Hometime","DistTravelled","RoG","MaxDiam","MaxHomeDist","SigLocsVisited","AvgFlightLen","StdFlightLen","AvgFlightDur","StdFlightDur","ProbPause","SigLocEntropy","MinsMissing","CircdnRtn","WkEndDayRtn")
  rownames(outmat)=daystr_v
  for(i in 1:length(daystr_v)){
    if(length(subsetinds_v[[i]])==0){next}
    submat=matrix(mobmat[subsetinds_v[[i]],],ncol=7)
    if(submat[1,1]==2 && submat[1,4]<subsetstarttime_v[[i]]){
      submat[1,4]=subsetstarttime_v[[i]]
    }
    if(submat[nrow(submat),1]==2 && submat[nrow(submat),7]>subsetstarttime_v[[i]]+60*60*24){
      submat[nrow(submat),7]=subsetstarttime_v[[i]]+60*60*24
    }
    submatmiss=matrix(mobmatmiss[subsetindsmiss_v[[i]],],ncol=7)
    if(submatmiss[1,1]==4 && submatmiss[1,4]<subsetstarttime_v[[i]]){
      submatmiss[1,4]=subsetstarttime_v[[i]]
    }
    if(submatmiss[nrow(submatmiss),1]==4 && submatmiss[nrow(submatmiss),7]>subsetstarttime_v[[i]]+60*60*24){
      submatmiss[nrow(submatmiss),7]=subsetstarttime_v[[i]]+60*60*24
    }
    if(nrow(submat)==0 || length(which(slout$home==1))==0){
      outmat[i,]=c(rep(NA,12),1440,rep(NA,2))
    }else{
      outmat[i,1]=Hometime(submat,slout,CENTERRAD=200)    
      outmat[i,2]=DistanceTravelled(submat)
      outmat[i,3]=RadiusOfGyration(submat,ITRVL)
      outmat[i,4]=MaxDiam(submat)
      outmat[i,5]=MaxHomeDist(submat,homex,homey)
      outmat[i,6]=SigLocsVisited(submat,slout,CENTERRAD)
      outmat[i,7]=AvgFlightLen(submat)
      outmat[i,8]=StdFlightLen(submat)
      outmat[i,9]=AvgFlightDur(submat)
      outmat[i,10]=StdFlightDur(submat)
      outmat[i,11]=ProbPause(submat)
      outmat[i,12]=SigLocEntropy(submat,slout,CENTERRAD)
      outmat[i,13]=MinsMissing(submatmiss)
      DRIout=DailyRoutineIndex(i,mobmat,subsetinds_v,subsetdayofweek_v,subsetstarttime_v,tz,CENTERRAD)
      outmat[i,14]=DRIout$cscore
      outmat[i,15]=DRIout$wscore      
    }
  }
  return(list(outmat,slout))
}
