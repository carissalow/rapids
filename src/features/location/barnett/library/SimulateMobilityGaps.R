SimulateMobilityGaps <-
function(suboutmat,obj,wtype="TL",spread_pars=c(1,10)){
  ind11=obj$ind11;ind12=obj$ind12;fd=obj$fd;ft=obj$ft;fts=obj$fts;fa=obj$fa;pt=obj$pt;pts=obj$pts;allts=obj$allts;phatall=obj$phatall;fxs=obj$fxs;fys=obj$fys;pxs=obj$pxs;pys=obj$pys;allxs=obj$allxs;allys=obj$allys
  if(nrow(suboutmat)==0){
    return(suboutmat)
  }
  foutmat=c()
  for(i in 1:nrow(suboutmat)){
    if(suboutmat[i,1]==1){
      curx=suboutmat[i,5]
      cury=suboutmat[i,6]
      foutmat=rbind(foutmat,suboutmat[i,])
    }else  if(suboutmat[i,1]<=3){
      curx=suboutmat[i,2] 
      cury=suboutmat[i,3]
      foutmat=rbind(foutmat,suboutmat[i,])
    }
    if(suboutmat[i,1]==4 && i>1 && i<nrow(suboutmat)){
      varmult=1
      while(TRUE){
        fw=dnorm((fts-mean(c(suboutmat[i,4],suboutmat[i,7])))/(varmult*(suboutmat[i,7]-suboutmat[i,4])))
        pw=dnorm((pts-mean(c(suboutmat[i,4],suboutmat[i,7])))/(varmult*(suboutmat[i,7]-suboutmat[i,4])))
        allw=dnorm((allts-mean(c(suboutmat[i,4],suboutmat[i,7])))/(varmult*(suboutmat[i,7]-suboutmat[i,4])))
        if(length(pts)>0 && length(fts)>0 && sum(fw)>0 && sum(pw)>0){break}
        if(length(pts)==0 && length(fts)>0 && sum(fw)>0){break}
        if(length(fts)==0 && length(pts)>0 && sum(pw)>0){break}
        varmult=varmult*2
      }
      s11=sum(allw[ind11],na.rm=T)
      s12=sum(allw[ind12],na.rm=T)
      if(s11+s12==0){phatcur=phatall}else{phatcur=s12/(s11+s12)}
      if(wtype=="LI"){
        foutmat=rbind(foutmat,c(1,curx,cury,suboutmat[i,4],suboutmat[i+1,2],suboutmat[i+1,3],suboutmat[i,7]))
      }else{
        rbout=matrix(RandomBridge(x0=curx,y0=cury,x1=suboutmat[i+1,2],y1=suboutmat[i+1,3],t0=suboutmat[i,4],t1=suboutmat[i,7],fd=fd,ft=ft,fts=fts,fa=fa,fw=fw,probp=phatcur,pt=pt,pts=pts,pw=pw,allts=allts,allw=allw,ind11=ind11,ind12=ind12,i_ind=i,pxs=pxs,pys=pys,fxs=fxs,fys=fys,allxs=allxs,allys=allys,wtype=wtype,canpause=suboutmat[i-1,1]==1,niter=100,spread_pars=spread_pars),ncol=7)
        foutmat=rbind(foutmat,rbout)        
      }
    }
  }  
  return(foutmat)
}
