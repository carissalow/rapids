MaxHomeDist <-
function(mat,homex,homey){
  IDmv=which(mat[,1]<=2)
  if(length(IDmv)==0){return(NA)}
  dfhome=rep(NA,length(IDmv))
  for(i in 1:length(IDmv)){
    dfhome[i]=sqrt((mat[IDmv[i],2]-homex)^2+(mat[IDmv[i],3]-homey)^2)
  }
  return(max(dfhome))
}
