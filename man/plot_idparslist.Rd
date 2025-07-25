% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/secsse_utils.R
\name{plot_idparslist}
\alias{plot_idparslist}
\title{function to visualize the structure of the idparslist}
\usage{
plot_idparslist(idparslist, state_names)
}
\arguments{
\item{idparslist}{idparslist setup list}

\item{state_names}{names of all states (including concealed states)}
}
\value{
list with two ggplot objects: plot_qmat which visualizes the
q_matrix, and plot_lambda which visualizes the lambda matrices.
}
\description{
function to visualize the structure of the idparslist
}
\examples{
idparslist <- list()
focal_matrix <-
secsse::create_default_lambda_transition_matrix(state_names = c("1", "2"),
                                                model = "CR")
idparslist[[1]] <- 
  secsse::create_lambda_list(state_names = c("1", "2"),
                             num_concealed_states = 2,
                             transition_matrix = focal_matrix,
                             model = "CR")
idparslist[[2]] <- secsse::create_mu_vector(state_names = c("1", "2"),
                                            num_concealed_states = 2,
                                            model = "CR",
                                            lambda_list = idparslist[[1]])
shift_mat <- secsse::create_default_shift_matrix(state_names = c("1", "2"),
                                                 num_concealed_states = 2,
                                                 mu_vector = idparslist[[2]])
idparslist[[3]] <- secsse::create_q_matrix(state_names = c("1", "2"),
                                           num_concealed_states = 2,
                                           shift_matrix = shift_mat,
                                           diff.conceal = FALSE)
p <- plot_idparslist(idparslist, state_names = names(idparslist[[1]]))
p$plot_lambda
p$plot_qmat
}
