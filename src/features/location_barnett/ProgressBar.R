ProgressBar <-
function (maxn, ind){
  if(maxn<51){
    if (ind == 1) {
      cat("|1%--------------------50%--------------------100%|\n")
      cat("|")
      return()
    }else{
      numprint=floor(50*ind/maxn)-floor(50*(ind-1)/maxn)
      if(numprint>0){
        for(i in 1:numprint){
          cat("|")
        }
      }
    }
  }else{
    if (ind == 1) {
      cat("|1%--------------------50%--------------------100%|\n")
      cat("|")
      return()
    }
    if (maxn == ind) {
      cat("|\n")
      return()
    }
    if (floor(50 * ind/maxn) != floor(50 * (ind - 1)/maxn)) {
      cat("|")
    }    
  }
}
