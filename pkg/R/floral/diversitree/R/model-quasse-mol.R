## TODO:
## 1. check that time is entirely completed.  This is important (done)
## 2. Remove saveOut() if possible for a 7% speedup, mostly by losing
##    the method dispatch on the transpose operation.
## 3. Add in the recursive underflow checking.  This might be hard to
##    do with the f.hi/f.lo combination, but there might be some hints
##    in the bisse.td code on how to do this nicely.
## 4. Tunable accuracy, probably via control.

make.branches.quasse.mol <- function(control) {
  nx <- control$nx
  dx <- control$dx
  r <- control$r
  tc <- control$tc
  atol <- control$atol
  rtol <- control$rtol
  
  f.hi <- make.pde.quasse.mol(nx*r, dx/r, 2L, atol, rtol)
  f.lo <- make.pde.quasse.mol(nx,   dx,   2L, atol, rtol)
  make.branches.quasse(f.hi, f.lo, control)
}

make.pde.quasse.mol <- function(nx, dx, nd, atol, rtol) {
  diffusion.scal <- 1/(2 * dx^2)

  function(y, len, pars, t0, method="lsodes", ...) {
    pars <- c(pars$lambda, pars$mu, pars$lambda + pars$mu,
              pars$drift, pars$diffusion * diffusion.scal)
    
    ans <- ode.1D(y, c(t0, t0+len), "derivs_quasse_mol", pars,
                  initfunc="initmod_quasse_mol",
                  nspec=nd, dimens=nx,
                  method=method, dllname="diversitree", atol=atol,
                  rtol=rtol)
    if ( abs(ans[length(len)+1,1] - (t0 + len)) > 1e-8 )
      stop("Integration error: integration stopped prematurely")
    
    matrix(ans[-1,-1], nx, nd)
  }
}
