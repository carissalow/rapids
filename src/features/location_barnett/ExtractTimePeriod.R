ExtractTimePeriod <-
function(tstart,tend,out){
  tstart = as.POSIXct(tstart)
  tend = as.POSIXct(tend)
  INDs=intersect(which(apply(out[,c(4,7)],1,function(x) max(x,na.rm=T))>=tstart),which(out[,4]<=tend))
  suboutmat=out[INDs,]  
  return(suboutmat)
}
