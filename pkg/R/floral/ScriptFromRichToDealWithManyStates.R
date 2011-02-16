#code from rich to do musse with many states


 ## With 7 binary traits, you have 2^7 = 128 states.  Simulate a tree
 ## so that we have something to make a likelihood function from.
 set.seed(1)
 k <- 128
 pars <- rep(c(.1, 0, .03), c(k, k, k*(k-1)))
 set.seed(1)
 phy <- tree.musse(pars, 100, x0=1)

 ## Here's the likelihood function with the collided argument names.
 lik <- make.musse(phy, phy$tip.state, k, strict=FALSE)
 argnames(lik) # 16512 arguments!

 ## Rather than work with 1..128, I'm guessing you have the traits
 ## scored as 0/1 for each -- this just makes binary labels like
 ## that.  (Here trait 1 cycles slowest, 7 fastest).  But this could
 ## be anything, including as.character(1:128)
 states <- apply(expand.grid(rep(list(0:1), 7))[7:1],
                 1, paste, collapse="")

 ## The only tricky part here is coming up with the correct names for
 ## the q matrix, as this should be by row and skipping diagonal
 ## elements, (q12, q13, ..., etc)
 tmp <- matrix(paste("q", rep(states, each=k), states, sep="."), k, k)
 diag(tmp) <- NA
 str.q <- tmp[!is.na(tmp)]
 str <- c(paste("lambda", states, sep="."),
          paste("mu", states, sep="."),
          str.q)

 ## The argument names of the likelihood function can then be
 ## directly replaced by this string.
 argnames(lik) <- str

 ## Shows the new names:
 argnames(lik)

 ## These will propagate through to constrain/find.mle/mcmc, etc.
 ## However, for optimisation/mcmc you should either use an unnamed
 ## starting vector, or name it with names(p) <- argnames(lik) to
 ## make the names stick.

 ## For constraining, one thing I've been making use of is generating
 ## constraint lists by pattern matching, as for functions with large
 ## numbers of parameters I make mistakes a lot.  For example,
 formulae.regexp <- function(pat, lik) {
   names <- argnames(lik)
   lhs <- as.character(pat)[[2]]
   rhs <- as.character(pat)[[3]]
   lhs.v <- setdiff(grep(lhs, argnames(lik), value=TRUE), rhs)
   lapply(sprintf("%s ~ %s", lhs.v, rhs), as.formula)
 }

 ## This uses regular expressions (not globs, but see ?glob2rx), but
 ## you can do things like
 formulae <- formulae.regexp("mu.*" ~ mu.0000000, lik)
 lik.constantmu <- constrain(lik, formulae=formulae)
 ## to set all the mu values to be the same across all states.

 ## Similarly, you could disallow all double moves (a la Pagel 1994)
 ## by doing:
 tmp <- expand.grid(rep(list(0:1), 7))[7:1]
 ok <- as.matrix(dist(tmp, "man")) == 1
 diag(ok) <- NA
 ## (surprisingly slow - it seems R's formulae objects are not the
 ## fastest things to work with in large numbers)
 formulae2 <- lapply(sprintf("%s ~ 0", str.q[!ok[!is.na(ok)]]),
                     as.formula)
 ## (sorry - also *really* slow, same reason)
 lik2 <- constrain(lik, formulae=c(formulae, formulae2))
 length(argnames(lik2)) # much smaller: 1025 arguments now...