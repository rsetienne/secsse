#' Print likelihood for initial parameters
#'
#' @inheritParams default_params_doc
#' @param initloglik A numeric with the value of loglikehood obtained prior to
#'   optimisation. Only used internally.
#'
#' @return Invisible `NULL`. Prints a `message()` to the console with the
#'   initial loglikelihood if `verbose >= 1`
#' @noRd
print_init_ll <- function(initloglik,
                          verbose) {
  if (isTRUE(verbose >= 1)) {
    init_ll_msg1 <- "Calculating the likelihood for the initial parameters."
    init_ll_msg2 <-
      paste0("The loglikelihood for the initial parameter values is ",
             initloglik)
    init_ll_msg3 <- c("Optimizing the likelihood - this may take a while.")
    message(paste(init_ll_msg1, init_ll_msg2, init_ll_msg3, sep = "\n"))
  }

  invisible(NULL)
}
