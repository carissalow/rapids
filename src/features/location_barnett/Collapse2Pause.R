Collapse2Pause <-
function(mat){
  cent=colMeans(mat[,2:3],na.rm=TRUE)
  if(!is.na(mat[nrow(mat),7])){
    return(c(2,cent[1],cent[2],mat[1,4],NA,NA,mat[nrow(mat),7]))
  }else{
    return(c(2,cent[1],cent[2],mat[1,4],NA,NA,mat[nrow(mat),4]))    
  }
}
