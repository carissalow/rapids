MaxRadius <-
function(mat){
  cent=colMeans(mat,na.rm=TRUE)
  return(max(apply(mat,1,function(x) sqrt((x[1]-cent[1])^2+(x[2]-cent[2])^2)),na.rm=TRUE))
}
