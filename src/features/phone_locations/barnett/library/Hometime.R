Hometime <-
function(mat,slout,CENTERRAD){
  IDhome=which(slout$home==1)
  xcenter=slout[IDhome,1]
  ycenter=slout[IDhome,2]
  tottime=0
  for(i in 1:nrow(mat)){
    if(mat[i,1]==2 && sqrt((mat[i,2]-xcenter)^2+(mat[i,3]-ycenter)^2)<CENTERRAD){ ### error missing value where true/false needed
      tottime=tottime+mat[i,7]-mat[i,4]
    }
  }
  return(tottime/60)
}
