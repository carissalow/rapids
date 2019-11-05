MaxDiam <-
function(mat){
  IDmv=which(mat[,1]<=2)
  if(length(IDmv)<2){return(0)}
  return(max(dist(mat[IDmv,2:3])))
}
