ybi <- function(Xn, Yn, hn, B0) {
  
  p = nrow(B0)
  d = ncol(B0)
  
  # --- get handle to objective function and derivative ......................
  tx = function(x) { matrix(x, p, d) }
  
  # --- get handle to objective function and derivative ......................
  Fhandle = F(F4DDS_binary_cpp, Xn, Yn, hn)
  dFhandle = dF(dF4DDS_binary_cpp, Xn, Yn, hn)
  mod = Module("ManifoldOptim_module", PACKAGE = "ManifoldOptim")
  prob = new(mod$RProblem, Fhandle, dFhandle)
  
  # B0 = ini_DDS(Xn,Yn,d);
  mani.defn = get.grassmann.defn(p, d)
  mani.params = get.manifold.params(IsCheckParams = TRUE)
  solver.params = get.solver.params(IsCheckParams = TRUE)
  x0 = as.numeric(B0)
  res = manifold.optim(prob, mani.defn, method = "RTRSR1",
                       mani.params = mani.params, solver.params = solver.params, x0 = x0)
  Bhat = tx(res$xopt)
  
  return(Bhat)
}
