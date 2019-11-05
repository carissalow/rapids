DayDist <-
function(i1,i2,mobmat,subsetinds_v,subsetstarttime_v,tz,CENTERRAD){
  chkpts=subsetstarttime_v[i1]+seq(from=30*60,by=60*60,length.out=24)
  chkpts2=subsetstarttime_v[i2]+seq(from=30*60,by=60*60,length.out=24)
  mat1=matrix(mobmat[subsetinds_v[[i1]],],ncol=7)
  mat2=matrix(mobmat[subsetinds_v[[i2]],],ncol=7)
  SamePlace = rep(NA,24)
  for(i in 1:24){
    loc1=LocationAt(mat1,chkpts[i])
    if(is.null(loc1)){
      next
    }
    IDs=c(which(abs(mat2[,4]-chkpts[i])%%(60*60*24)<30*60),which(abs(mat2[,4]-chkpts[i])%%(60*60*24)>(60*60*24)-30*60))
    if(length(IDs)>0){
      CanBe0=FALSE
      for(j in 1:length(IDs)){
        if(mat2[IDs[j],1]<=3){
          CanBe0 = TRUE
          if(sqrt((mat2[IDs[j],2]-loc1$x)^2+(mat2[IDs[j],3]-loc1$y)^2)<CENTERRAD){
            SamePlace[i]=1
          }
        }
      }
      if(CanBe0 && is.na(SamePlace[i])){
        SamePlace[i]=0
      }
    }else{
      for(j in 1:nrow(mat2)){
        if(!is.na(mat2[j,4])&&!is.na(mat2[j,7])&& mat2[j,4]<chkpts2[i] && mat2[j,7]>chkpts2[i]){
          if(mat2[j,1]!=4 && sqrt((mat2[j,2]-loc1$x)^2+(mat2[j,3]-loc1$y)^2)<CENTERRAD){
            SamePlace[i]=1
          }else{
            SamePlace[i]=0
          }
          break
        }
      }
    }
  }
  if(length(which(!is.na(SamePlace)))==0){return(NA)}
  return(mean(SamePlace,na.rm=T))
}
