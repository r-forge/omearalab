\name{make.prior}
\alias{make.prior}
\alias{make.prior.exponential}

\title{Simple Prior Functions}

\description{Functions for generating prior functions for use with
  \code{\link{mcmc}}, etc.}

\usage{
make.prior.exponential(r)
}

\arguments{
  \item{r}{Scalar or vector of rate parameters}
}

\details{
  The prior probability distribution has probability density
  \deqn{\sum_i r_i e^{-r_i x_i}}{sum(r[i]*exp(-r[i]*x[i]))} where the
  \eqn{i} denotes the \eqn{i}th parameter.  If \code{r} is a scalar,
  then the same rate is used for all parameters.

  This function returns a function that may be used as the \code{prior}
  argument to the likelihood functions returned elsewhere.
}

\keyword{models}
\author{Richard G. FitzJohn}
