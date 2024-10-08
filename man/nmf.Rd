% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hnmf.R
\name{nmf}
\alias{nmf}
\title{NMF Iteration}
\usage{
nmf(
  catalog,
  sig = NULL,
  expos = NULL,
  denovo = TRUE,
  K = NULL,
  nrun = 10,
  verbose = 1,
  progress.bar = FALSE,
  Itmax = 1e+05,
  Tol = 1e-05,
  a = 10,
  nprint = 100,
  useC = FALSE,
  alpha = 1,
  initializer = "random",
  fix.h = NA,
  clustering = FALSE,
  logLik.cut = 0.05
)
}
\arguments{
\item{catalog}{Catalog matrix}

\item{denovo}{Full iteration of both signature and exposure; if \code{FALSE} (\code{hnmf}),
only the exposure is fit (requires \code{signat} input).}

\item{K}{Number of signatures.}

\item{nrun}{Number of independent runs to generate}

\item{verbose}{Verbosity level}

\item{progress.bar}{Display progress bar.}

\item{Itmax}{Maximum no. of iteration.}

\item{Tol}{Tolerance for checking convergence.}

\item{fix.h}{Binary flags for exposure to keep fixed during inference (for deviance estimate)}

\item{clustering}{Perform clustering of local maxima}

\item{logLik.cut}{Log likelihood fractional cutoff with respect to global maximum for clustering}

\item{signat}{Signature matrix}
}
\value{
Object of class \code{DeepSig}.
}
\description{
Non-negative matrix factorization with either fixed or variable signature matrix.
}
\details{
If \code{denovo = TRUE}, the full NMF is performed with both \code{W} and \code{H} matrix 
determined de novo. If \code{denovo = FALSE}, \code{W} is fixed as input and only \code{H} 
is determined.
}
