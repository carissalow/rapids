GuessPause <-
function(mat,mindur=300,r=75){
  cat("Inferring pauses...\n")
  flatmat=c()
  collapse=FALSE
  inds=1
  incr=1
  tcur=mat[inds,4]
  while(TRUE){
    if(!is.na(mat[inds+incr,7])){
      tnex=mat[inds+incr,7]
    }else{
      tnex=mat[inds+incr,4]
    }
    if(tnex-tcur>=mindur){
      if(mat[inds+incr,1]==1){
        maxr=MaxRadius(rbind(mat[inds:(inds+incr),2:3],mat[inds:(inds+incr),5:6]))
      }else{
        maxr=MaxRadius(mat[inds:(inds+incr),2:3])        
      }
      if(maxr>r && !collapse){
        inds=inds+1
        tcur=mat[inds,1]
        incr=0
      }
      if(maxr>r && collapse){
        if(mat[inds+incr-1,1]==4){
          if(inds<inds+incr-2){
            flatmat=rbind(flatmat,c(inds,inds+incr-2))                                
          }
        }else{
          flatmat=rbind(flatmat,c(inds,inds+incr-1))          
        }
        inds=inds+incr
        tcur=mat[inds,1]
        incr=0
        collapse=FALSE
      }
      if(maxr<=r){
        collapse=TRUE
      }
    }
    incr=incr+1
    if(inds+incr>nrow(mat)){
      if(maxr<=r && tnex-tcur>=mindur){
        flatmat=rbind(flatmat,c(inds,inds+incr-1))
      }
      break
    }
  }
  if(nrow(flatmat)==0){
    return(mat)
  }else{
    outmat=c()
    if(flatmat[1,1]>1){
      outmat=mat[1:(flatmat[1,1]-1),]
    }
    for(i in 1:nrow(flatmat)){
      #ProgressBar(nrow(flatmat),i)
      outmat=rbind(outmat,Collapse2Pause(mat[flatmat[i,1]:flatmat[i,2],]))
      if(i<nrow(flatmat) && flatmat[i,2]<flatmat[i+1,1]-1){
        outmat=rbind(outmat,mat[(flatmat[i,2]+1):(flatmat[i+1,1]-1),])
      }
    }    
    if(flatmat[nrow(flatmat),2]<nrow(mat)){
      outmat=rbind(outmat,mat[(flatmat[nrow(flatmat),2]+1):nrow(mat),])
    }
  }
  if(nrow(outmat)==1){
    rownames(outmat)=NULL
    colnames(outmat)=c("Code","x0","y0","t0","x1","y1","t1")
    return(outmat)
  }
  # Group together adjacent pauses, averaging their location
  flatmat2=c()
  outmat2=c()
  collapse=FALSE
  for(i in 2:nrow(outmat)){
    if(outmat[i,1]!=2 && !collapse){
      next
    }else if(outmat[i,1]!=2 && collapse){
      collapse=FALSE
      flatmat2=rbind(flatmat2,c(cstart,i-1))
    }else  if(outmat[i,1]==2 && outmat[i-1,1]==2 && !collapse){
      cstart=i-1
      collapse=TRUE
    }else if(outmat[i,1]==2 && collapse){
      next
    }
  }
  if(collapse && outmat[nrow(outmat),1]==2){
    flatmat2=rbind(flatmat2,c(cstart,nrow(outmat)))
  }
  if(is.null(flatmat2)){
    outmat2=outmat
  }else{
    flatmat2=matrix(flatmat2,ncol=2)
    outmat2=c()
    if(flatmat2[1,1]>1){
      outmat2=outmat[1:(flatmat2[1,1]-1),]
    }
    for(i in 1:nrow(flatmat2)){
      outmat2=rbind(outmat2,Collapse2Pause(outmat[flatmat2[i,1]:flatmat2[i,2],]))
      if(i<nrow(flatmat2) && flatmat2[i,2] < flatmat2[i+1,1]-1){
        outmat2=rbind(outmat2,outmat[(flatmat2[i,2]+1):(flatmat2[i+1,1]-1),])
      }
    }
    if(flatmat2[nrow(flatmat2),2]<nrow(outmat)){
      outmat2=rbind(outmat2,outmat[(flatmat2[nrow(flatmat2),2]+1):nrow(outmat),])
    }
  }
  # Set flight endpoints equal to pause endpoints
  if(outmat2[1,1]==1 && outmat2[2,1]==2){
    outmat2[1,5]=outmat2[2,2]
    outmat2[1,6]=outmat2[2,3]
  }
  for(i in 2:(nrow(outmat2)-1)){
    if(outmat2[i,1]==1 && outmat2[i-1,1]==2){
      outmat2[i,2]=outmat2[i-1,2]
      outmat2[i,3]=outmat2[i-1,3]
    }
    if(outmat2[i,1]==1 && outmat2[i+1,1]==2){
      outmat2[i,5]=outmat2[i+1,2]
      outmat2[i,6]=outmat2[i+1,3]
    }
  }
  if(outmat2[nrow(outmat2)-1,1]==2 && outmat2[nrow(outmat2),1]==1){
    outmat2[nrow(outmat2),2]=outmat2[nrow(outmat2)-1,2]
    outmat2[nrow(outmat2),3]=outmat2[nrow(outmat2)-1,3]
  }
  rownames(outmat2)=NULL
  colnames(outmat2)=c("Code","x0","y0","t0","x1","y1","t1")
  return(outmat2)
}
