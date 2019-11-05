SigLocsVisited <-
function(mat,slout,CENTERRAD){
  places_visited = rep(0,nrow(slout))
  for(i in 1:nrow(mat)){
    if(mat[i,1]<=3){
      for(j in 1:nrow(slout)){
        if(sqrt((slout[j,1]-mat[i,2])^2 + (slout[j,2]-mat[i,3])^2)<CENTERRAD){
          places_visited[j]=1
        }
      }
    }
  }
  return(sum(places_visited))
}
