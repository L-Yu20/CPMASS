dF <- function(fun_handle,Xn,Yn,hn){
  p <- ncol(Xn)
  
  df <- function(W){
    W <- matrix(W, nrow = p)
    diffval <- fun_handle(W,Xn,Yn,hn)
    return(diffval)
  }
  
  return(df)
}

