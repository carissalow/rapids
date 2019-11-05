SigLocEntropy <-
function(mat,slout,CENTERRAD){
  tp = rep(0,nrow(slout))
  for(i in 1:nrow(mat)){
    if(mat[i,1]==2){
      for(j in 1:nrow(slout)){
        if(sqrt((slout[j,1]-mat[i,2])^2+(slout[j,2]-mat[i,3])^2)<CENTERRAD){
          tp[j] = tp[j] + mat[i,7]-mat[i,4]
        }
      } 
    }
  }
  tot=0
  if(sum(tp)==0){return(0)}
  for(i in 1:nrow(slout)){
    p=tp[i]/sum(tp)
    if(p>0){
      tot=tot-p*log(p)
    }
  }
  return(tot)
}
