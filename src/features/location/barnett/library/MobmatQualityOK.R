MobmatQualityOK <-
function(mobmat,obj){
  msg=""
  if(!is.matrix(mobmat)){
    msg = "Mobmat not a matrix. Removing individual from analysis.\n"
  }
  if(length(obj$ID1)==0){
    msg= "No flights in mobmat. Removing individual from analysis.\n"
  }
  if(length(obj$ID2)==0){
    msg= "No pauses in mobmat. Removing individual from analysis.\n"
  }
  return(msg)
}
