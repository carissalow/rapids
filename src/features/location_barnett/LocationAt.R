LocationAt <-
function(mat,tt){
  for(i in 1:nrow(mat)){
    if(mat[i,1]<=2){
      if(mat[i,4]<=tt && mat[i,7]>=tt){
        return(list('x'=mat[i,2],'y'=mat[i,3]))
      }
    }else if(mat[i,1]==3){
      if(mat[i,4]==tt){
        return(list('x'=mat[i,2],'y'=mat[i,3]))
      }      
    }
  }
  return(NULL)
}
