AvgFlightLen <-
function(mat){
  num=0
  tot=0
  for(i in 1:nrow(mat)){
    if(mat[i,1]==1){
      tot=tot+sqrt((mat[i,5]-mat[i,2])^2+(mat[i,6]-mat[i,3])^2)
      num=num+1
    }
  }
  if(num==0){return(0)}
  return(tot/num)
}
