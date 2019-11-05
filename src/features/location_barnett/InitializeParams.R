InitializeParams <-
function(out){
  ID1=which(out[,1]==1)
  ID2=which(out[,1]==2)
  ID3=which(out[,1]==3)
  ID4=which(out[,1]==4)
  
  # probability of a pause after a flight
  ID1p1=ID1+1    
  if(length(ID1)>0 && ID1[length(ID1)]==nrow(out)){  
    ID1p1=ID1p1[-length(ID1p1)]
  }
  allts=apply(out,1,function(xx) mean(xx[c(4,7)]))
  allxs=out[,2]
  allys=out[,3]
  ind11=ID1p1[which(out[ID1p1,1]==1)]
  ind12=ID1p1[which(out[ID1p1,1]==2)]
  l1=length(ind11)
  l2=length(ind12)
  if(l1+l2>0){
    phatall=l2/(l1+l2)    
  }
  if(l1+l2==0){phatall=length(ID2)/(length(ID1)+length(ID2))}
  #flight distances
  fd=apply(out[ID1,],1,function(xx) sqrt((xx[2]-xx[5])^2+(xx[3]-xx[6])^2))
  
  # flight times: ft
  ft=apply(out[ID1,],1,function(xx) (xx[7]-xx[4]))
  fxs=out[ID1,2]
  fys=out[ID1,3]
  # flight angles range [0,2pi]: fa
  #fa=apply(out[ID1,],1,function(xx) atan((xx[6]-xx[3])/(xx[5]-xx[2]))-((sign(xx[6]-xx[3])-1)/2)*pi)
  fa=rep(0,length(ID1))
  yvals=out[ID1,6]-out[ID1,3]
  xvals=out[ID1,5]-out[ID1,2]
  IDyg0=which(yvals>=0)
  IDxg0=which(xvals>=0)
  IDyl0=which(yvals<0)
  IDxl0=which(xvals<0)  
  IDgg=intersect(IDyg0,IDxg0)
  IDlg=intersect(IDyg0,IDxl0)
  IDgl=intersect(IDyl0,IDxg0)
  IDll=intersect(IDyl0,IDxl0)
  fa[IDgg]=atan(yvals[IDgg]/xvals[IDgg])
  fa[IDgl]=atan(yvals[IDgl]/xvals[IDgl])+2*pi
  fa[IDlg]=atan(yvals[IDlg]/xvals[IDlg])+pi
  fa[IDll]=atan(yvals[IDll]/xvals[IDll])+pi
  # flight time stamps: fts
  fts=out[ID1,4]
  
  # pause times
  pt=apply(matrix(out[ID2,],ncol=7),1,function(xx) xx[7]-xx[4])
  pxs=out[ID2,2]
  pys=out[ID2,3]
  #pause time stamp: pts
  pts=out[ID2,4]  
  return(list(ID1=ID1,ID2=ID2,ID3=ID3,ID4=ID4,ID1p1=ID1p1,allts=allts,ind11=ind11,ind12=ind12,phatall=phatall,fd=fd,ft=ft,fa=fa,fts=fts,pt=pt,pts=pts,fxs=fxs,fys=fys,pxs=pxs,pys=pys,allxs=allxs,allys=allys))  
}
