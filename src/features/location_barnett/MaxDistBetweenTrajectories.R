MaxDistBetweenTrajectories <-
function(mat1,mat2,t_gap=1){
  mat1=matrix(mat1,ncol=7);mat2=matrix(mat2,ncol=7)
  t0=mat1[1,4];t1=mat2[nrow(mat2),7]
  t_mesh = seq(t0,t1,t_gap)
  d_v = rep(0,length(t_mesh))
  for(i in 1:length(t_mesh)){
    ID1=intersect(which(mat1[,4]<=t_mesh[i]),which(mat1[,7]>=t_mesh[i]))[1]
    if(mat1[ID1,1]==2){
      x1=mat1[ID1,2]
      y1=mat1[ID1,3]
    }else{
      w1=(t_mesh[i]-mat1[ID1,4])/(mat1[ID1,7]-mat1[ID1,4])
      x1=mat1[ID1,2]*(1-w1)+mat1[ID1,5]*w1
      y1=mat1[ID1,3]*(1-w1)+mat1[ID1,6]*w1
    }
    ID2=intersect(which(mat2[,4]<=t_mesh[i]),which(mat2[,7]>=t_mesh[i]))[1]
    if(mat2[ID2,1]==2){
      x2=mat2[ID2,2]
      y2=mat2[ID2,3]
    }else{
      w2=(t_mesh[i]-mat2[ID2,4])/(mat2[ID2,7]-mat2[ID2,4])
      x2=mat2[ID2,2]*(1-w2)+mat2[ID2,5]*w2
      y2=mat2[ID2,3]*(1-w2)+mat2[ID2,6]*w2
    }
    d_v[i] = sqrt((x1-x2)^2+(y1-y2)^2)
  }
  return(max(d_v))
}
