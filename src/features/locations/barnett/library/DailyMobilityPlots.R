DailyMobilityPlots <-
function(mobmat,obj,tz,filename){
  curdate=strsplit(as.character(as.POSIXct(mobmat[1,4],tz=tz,origin="1970-01-01"))," ")[[1]][1]
  curtime=strsplit(strsplit(as.character(as.POSIXct(mobmat[1,4],tz=tz,origin="1970-01-01"))," ")[[1]][2],":")
  subsetinds_v = list()
  daystr_v = c(curdate)
  dayind=1
  subsetinds = c(1)
  for(i in 2:nrow(mobmat)){
    nexdate=strsplit(as.character(as.POSIXct(mobmat[i,4],tz=tz,origin="1970-01-01"))," ")[[1]][1]
    if(nexdate == curdate && i<nrow(mobmat)){
      subsetinds=c(subsetinds,i)
    }else{
      subsetinds_v[[dayind]]=subsetinds
      dayind=dayind+1
      if(mobmat[i-1,1]==2 && (mobmat[i-1,7]-mobmat[i-1,4])/(60*60*24)>1){
        ii=1
        middate = strsplit(as.character(as.POSIXct(mobmat[i-1,4]+(60*60*24)*ii,tz=tz,origin="1970-01-01"))," ")[[1]][1]
        while(middate!=nexdate){
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
      if(curdate!=daystr_v[length(daystr_v)] && length(daystr_v)<length(subsetinds_v)){
        daystr_v=c(daystr_v,curdate)        
      }
    }
  }  
  plot.flights(mobmat,diminch=4,outfile=paste("FlightsPlot_full_",filename,".pdf",sep=""),xrang=plotlimits(mobmat)$xrang,yrang=plotlimits(mobmat)$yrang)
  for(i in 1:length(daystr_v)){
    submat=matrix(mobmat[subsetinds_v[[i]],],ncol=7)
    plot.flights(submat,diminch=4,outfile=paste("FlightsPlot_",daystr_v[i],"_ZOOMOUT_",filename,".pdf",sep=""),xrang=plotlimits(mobmat)$xrang,yrang=plotlimits(mobmat)$yrang)      
    plot.flights(submat,diminch=4,outfile=paste("FlightsPlot_",daystr_v[i],"_ZOOMIN_",filename,".pdf",sep=""))      
  }
}
