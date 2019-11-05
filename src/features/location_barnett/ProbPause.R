ProbPause <-
function(mat){
  tpause = 0 
  tflight = 0
  for(i in 1:nrow(mat)){
    if(mat[i,1]==1){
      tflight = tflight + mat[i,7]-mat[i,4]
    }
    if(mat[i,1]==2){
      tpause = tpause + mat[i,7]-mat[i,4]
    }
  }
  return(tpause/(tpause+tflight))
}
