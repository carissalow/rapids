SigLocs <-
function(mobmat,obj,CENTERRAD=125,MINPAUSETIME=600,tz=""){
  if(length(obj$ID2)==0){
    warning("No pauses in mobmat within function SigLocs!")
    return(NULL)
  }else if(length(obj$ID2)==1){
    outmat=data.frame('x'=mobmat[obj$ID2[1],2],'y'=mobmat[obj$ID2[1],3],'timepresent'=c(0),'home'=c(1))
    nrowfc=1
  }else{  
    ptred=floor(obj$pt/MINPAUSETIME)
    if(length(which(ptred>0))<2){
      warning("No pauses long enough in mobmat within function SigLocs!")
      return(NULL)
    }
    pmat=c()
    for(i in 1:length(obj$ID2)){
      if(ptred[i]>0){
        pmat=rbind(pmat,matrix(rep(mobmat[obj$ID2[i],2:3],ptred[i]),ncol=2,byrow=T))
      }
    }  
    kmeansk_v=2:length(which(ptred>0))
    lsfit = list()
    for(i in 1:length(kmeansk_v)){
      kmeansk = kmeansk_v[i]
      fit = kmeans(pmat,centers=kmeansk)
      lsfit[[i]]=fit
      if(min(dist(fit$centers))<CENTERRAD*2 || i==length(kmeansk_v)){
        if(i>1){
          kmeansk=kmeansk_v[i-1]
          fit = lsfit[[i-1]]      
        }
        break
      }
    }
    nrowfc = nrow(fit$centers)
    outmat=data.frame('x'=fit$centers[,1],'y'=fit$centers[,2],'timepresent'=rep(0,nrow(fit$centers)),'home'=rep(0,nrow(fit$centers)))
  }
  #Determine time spent at these significant locations
  for(i in 1:length(obj$ID2)){
    for(j in 1:nrowfc){
      if(sqrt((mobmat[obj$ID2[i],2]-outmat$x[j])^2+(mobmat[obj$ID2[i],3]-outmat$y[j])^2)<CENTERRAD){
        outmat[j,3]=outmat[j,3]+obj$pt[i]
      }
    }
  }
  #Determine which is home (where is the night spent)
  for(i in 1:length(obj$ID2)){
    xx=as.POSIXct((mobmat[obj$ID2[i],7]+mobmat[obj$ID2[i],4])/2,tz=tz,origin="1970-01-01")
    hourofday = as.numeric(strsplit(strsplit(as.character(xx),":")[[1]][1]," ")[[1]][2])
    if(hourofday>=21 || hourofday<6){
      for(j in 1:nrowfc){
        if(sqrt((mobmat[obj$ID2[i],2]-outmat$x[j])^2+(mobmat[obj$ID2[i],3]-outmat$y[j])^2)<CENTERRAD){
          outmat[j,4]=outmat[j,4]+obj$pt[i]
        }
      }      
    }
  }
  IDmax =order(outmat[,4],decreasing=T)[1] 
  outmat[,4]=rep(0,nrow(outmat))
  outmat[IDmax,4]=1
  outmat=outmat[order(outmat[,3],decreasing=T),]
  IDrm=which(outmat[,3]==0)
  if(length(IDrm)>0){
    outmat=outmat[-IDrm,]
  }
  rownames(outmat)=NULL
  return(outmat)
}
