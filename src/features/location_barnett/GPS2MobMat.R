GPS2MobMat <-
function(locations_df,itrvl=10,accuracylim=51,r=NULL,w=NULL,tint_m=NULL,tint_k=NULL){
  if(is.null(r)){
    r=sqrt(itrvl)
  }
  cat("Read GPS coordinates...\n")
  mat=as.matrix(locations_df)
  colnames(mat)=c("timestamp","latitude","longitude","altitude","accuracy")
  mat=data.frame(mat)
  mat = mat[order(mat[,1]),]
  mat=mat[which(mat$accuracy<accuracylim),]
  if(!is.null(tint_k) && !is.null(tint_m)){
    t0 = mat$timestamp[1]/1000;
    mat=mat[which((mat$timestamp/1000-t0)%%(tint_k+tint_m)<tint_k),]    
  }
  if(is.null(w)){
    w=mean(mat$accuracy)+itrvl
  }
  tstart=mat[1,1]/1000
  tend=mat[nrow(mat),1]/1000
  avgmat = matrix(NA,nrow=ceiling((tend-tstart)/itrvl)+2,ncol=4)
  IDam = 1
  count = 0
  nextline=c(1,tstart+itrvl/2,mat[1,2],mat[1,3]) 
  numitrvl=1
  cat("Collapse data within",itrvl,"second intervals...\n")
  for(i in 2:nrow(mat)){
    ProgressBar(nrow(mat)-1,i-1)
    if(mat[i,1]/1000<tstart+itrvl){
      nextline[3]=nextline[3]+mat[i,2]
      nextline[4]=nextline[4]+mat[i,3]
      numitrvl=numitrvl+1
    }else{
      nextline[3]=nextline[3]/numitrvl
      nextline[4]=nextline[4]/numitrvl
      #avgmat=rbind(avgmat,nextline)
      avgmat[IDam,]=nextline
      count=count+1
      IDam=IDam+1
      nummiss=floor((mat[i,1]/1000-(tstart+itrvl))/itrvl)
      if(nummiss>0){
        #avgmat = rbind(avgmat,c(4,tstart+itrvl/2,tstart+itrvl*(nummiss+1)+itrvl/2,NA))
        avgmat[IDam,] = c(4,tstart+itrvl/2,tstart+itrvl*(nummiss+1)+itrvl/2,NA)
        count=count+1
        IDam=IDam+1
      }
      tstart=tstart+itrvl*(nummiss+1)
      nextline[1]=1
      nextline[2]=tstart+itrvl/2
      nextline[3]=mat[i,2]
      nextline[4]=mat[i,3]
      numitrvl=1
    }
  }
  avgmat = avgmat[1:count,]
  avgmat=cbind(avgmat[,1:4],NA,NA)
  ID1 = which(avgmat[,1]==1)
  cat("Convert from Lat/Lon to X/Y...\n")
  obj=LatLong2XY(avgmat[ID1,3],avgmat[ID1,4])
  avgmat[ID1,5:6]=cbind(obj$x_v,obj$y_v)
  outmat=c()
  curind=1
  cat("Convert from X/Y to flights/pauses...\n")
  for(i in 1:nrow(avgmat)){
    ProgressBar(nrow(avgmat),i)
    if(avgmat[i,1]==4){
      outmat=rbind(outmat,ExtractFlights(avgmat[curind:(i-1),c(5,6,2)],r,w),
                   c(avgmat[i,1],NA,NA,avgmat[i,2],NA,NA,avgmat[i,3]))
      curind=i+1
    }
  }  
  if(curind<=nrow(avgmat)){  
    outmat=rbind(outmat,ExtractFlights(avgmat[curind:nrow(avgmat),c(4,3,2)],r,w))
  }
  rownames(outmat)=NULL
  colnames(outmat)=c("Code","x0","y0","t0","x1","y1","t1")
  return(outmat)
}
