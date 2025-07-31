F <- function(fun_handle,Xn,Yn,hn){
  
  # Generic function to compute the objective function F. The function first
  # sets a handle to the specific model function and fixes the parameters from
  # the sample needed for its computation. The handle fixed with those
  # parameters is then evaluated at a given value for argument W.
  # =========================================================================
  p <- ncol(Xn)
  
  f <- function(W){
    W <- matrix(W, nrow = p)
    fval <- fun_handle(W,Xn,Yn,hn)
    return(fval)
  }
  
  return(f)
}

