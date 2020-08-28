StdFlightLen <-
function(mat){
  ID1=which(mat[,1]==1)
  if(length(ID1)<=1){return(0)}
  try1=try(sd(as.numeric(sqrt((mat[ID1,6]-mat[ID1,3])^2+(mat[ID1,5]-mat[ID1,2])^2))),silent=TRUE)
  if(class(try1) == "try-error" || is.na(try1)){
    return(0)
  }else{
    return(try1)
  }
}
