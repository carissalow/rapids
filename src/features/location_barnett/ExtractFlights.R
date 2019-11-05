ExtractFlights <-
function(mat,r,w){
  out = c()  
  if(!is.matrix(mat)){
    out=matrix(c(3,mat[1],mat[2],mat[3],NA,NA,NA),nrow=1)
    colnames(out)=c("code","lon1","lat1","t1","lon2","lat2","t2")
    return(out)
  }
  nextline = rep(NA,6)
  nextline[1:3]=mat[1,]
  curind = 1
  while(TRUE){
    nexind = curind+1
    if(nexind==nrow(mat)){
      nextline[4:6]=mat[nexind,]
      out=rbind(out,nextline)
      break
    }
    while(TRUE){
      if(!IsFlight(mat[curind:nexind,],r,w)){
        break
      }
      nexind=nexind+1
      if(nexind>nrow(mat)){
        break
      }
    }
    if(nexind==curind+1 && curind!=nrow(mat)){
      nextline[3]=mat[nexind,3]
      mat=mat[-nexind,]
    }else{
      nextline[4:6]=mat[nexind-1,]
      out=rbind(out,nextline)
      nextline=rep(NA,6)
      curind=nexind-1
      nextline[1:3]=mat[curind,]
    }
    if(nexind>nrow(mat)){break}
  }
  outp=c()
  if(out[1,3]>min(mat[,3])){
    outp=rbind(outp,c(2,out[1,1],out[1,2],min(mat[,3]),NA,NA,out[1,3]))
    if(nrow(out)==1){
      outp[1,7]=out[1,6]
      return(outp)
    }
  }
  if(nrow(out)==1){
    outp=rbind(outp,c(1,out))
    return(outp)
  }
  for(i in 1:(nrow(out)-1)){
    outp=rbind(outp,c(1,out[i,]))
    if(out[i,6]<out[i+1,3]){
      outp=rbind(outp,c(2,out[i,4],out[i,5],out[i,6],NA,NA,out[i+1,3]))
    }
  }
  # convert flights with distance 0 into pauses
  IDf=which(outp[,1]==1)
  IDp=which((outp[IDf,2]-outp[IDf,5])^2+(outp[IDf,3]-outp[IDf,6])^2==0)
  if(length(IDp)>0){
    for(i in 1:length(IDp)){
      outp[IDf[IDp[i]],1]=2
      outp[IDf[IDp[i]],5]=NA
      outp[IDf[IDp[i]],6]=NA
    }
  }
  outp=rbind(outp,c(1,out[nrow(out),]))
  colnames(outp)=c("code","lon1","lat1","t1","lon2","lat2","t2")
  return(outp)
}
