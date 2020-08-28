DistanceTravelled <-
function(mat){
  dt=0
  for(i in 1:nrow(mat)){
    if(mat[i,1]==1){
      dt=dt+sqrt((mat[i,5]-mat[i,2])^2+(mat[i,6]-mat[i,3])^2)
    }
  }
  return(dt)
}
