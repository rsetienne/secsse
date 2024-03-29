% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/secsse_prep.R
\name{create_default_lambda_transition_matrix}
\alias{create_default_lambda_transition_matrix}
\title{Helper function to create a default lambda list}
\usage{
create_default_lambda_transition_matrix(
  state_names = c("0", "1"),
  model = "ETD"
)
}
\arguments{
\item{state_names}{vector of names of all observed states.}

\item{model}{used model, choice of \code{"ETD"} (Examined Traits Diversification),
\code{"CTD"} (Concealed Traits Diversification) or \code{"CR"} (Constant Rate).}
}
\description{
This function generates a generic lambda list, assuming no transitions
between states, e.g. a species of observed state 0 generates daughter
species with state 0 as well.
}
\examples{
lambda_matrix <-
     create_default_lambda_transition_matrix(state_names = c(0, 1),
                                             model = "ETD")
lambda_list <- create_lambda_list(state_names = c(0, 1),
                                  num_concealed_states = 2,
                                  transition_matrix = lambda_matrix,
                                  model = "ETD")
}
