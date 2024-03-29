% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/secsse_prep.R
\name{fill_in}
\alias{fill_in}
\title{Helper function to enter parameter value on their right place}
\usage{
fill_in(object, params)
}
\arguments{
\item{object}{lambda matrices, \code{q_matrix} or mu vector.}

\item{params}{parameters in order, where each value reflects the value
of the parameter at that position, e.g. \code{c(0.3, 0.2, 0.1)} will fill out
the value 0.3 for the parameter with rate identifier 1, 0.2 for the
parameter with rate identifier 2 and 0.1 for the parameter with rate
identifier 3.}
}
\value{
lambda matrices, \code{q_matrix} or mu vector with the correct values in
their right place.
}
\description{
Helper function to enter parameter value on their right place
}
