StdFlightDur <-
function(mat){
  ID1=which(mat[,1]==1)
  if(length(ID1)==0){return(0)}
  try1=try(sd(as.numeric(mat[ID1,7]-mat[ID1,4])),silent=TRUE)
  if(class(try1) == "try-error" || is.na(try1)){
    return(0)
  }else{
    return(try1)
  }
}
