MinsMissing <-
function(mat){
  tot=0
  for(i in 1:nrow(mat)){
    if(mat[i,1]==4){
      tot = tot+mat[i,7]-mat[i,4]
    }
  }
  return(tot/60)
}
