plotlimits <-
function(mat,defaultdist=100){
  xrang=range(c(mat[which(mat[,1]<=2),2],mat[which(mat[,1]<=1),5]))
  yrang=range(c(mat[which(mat[,1]<=2),3],mat[which(mat[,1]<=1),6]))
  if(xrang[2]==xrang[1]){
    xrang[2]=xrang[2]+defaultdist/2
    xrang[1]=xrang[1]-defaultdist/2
  }
  if(yrang[2]==yrang[1]){
    yrang[2]=yrang[2]+defaultdist/2
    yrang[1]=yrang[1]-defaultdist/2
  }
  return(list(xrang=xrang,yrang=yrang))
}
