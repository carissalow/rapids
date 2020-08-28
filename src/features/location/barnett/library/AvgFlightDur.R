AvgFlightDur <-
function(mat){
  num=0
  tot=0
  for(i in 1:nrow(mat)){
    if(mat[i,1]==1){
      tot=tot+mat[i,7]-mat[i,4]
      num=num+1
    }
  }
  if(num==0){return(0)}
  return(tot/num)
}
