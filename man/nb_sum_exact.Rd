% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nb_sum_exact.R
\name{nb_sum_exact}
\alias{nb_sum_exact}
\title{Implements Furman's exact PMF for the evaluation of the sum of arbitrary NB random variables. Called by other functions. Not intended to be run alone.}
\usage{
nb_sum_exact(phis, ps, n.terms = 1000, counts, n.cores = 1, tolerance = 0.001)
}
\arguments{
\item{phis}{Vector of individual dispersion parameters. Equivalent to 'size' in dnbinom.}

\item{ps}{Vector of individual probabilities of success.}

\item{n.terms}{The number of terms to include in the series for evaluating the PMF at a given number of counts. Defaults to 1000.}

\item{counts}{The vector of counts over which the PMF is evaluated.}

\item{n.cores}{The number of CPU cores to use in the evaluation. Allows parallelization.}

\item{tolerance}{The acceptable difference between the sum of the K distribution and 1.}
}
\value{
A numeric vector of probability densities.
}
\description{
Implements Furman's exact PMF for the evaluation of the sum of arbitrary NB random variables. Called by other functions. Not intended to be run alone.
}
\examples{
nb_sum_exact(ps = c(0.05, 0.44), phis = c(5, 8), counts = 0:500)

}
