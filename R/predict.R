#' @export
predict.ipriorBVS <- function(object, truth, ...) {
  pips <- get_pips(object)
  brier <- mean((pips - truth) ^ 2)
  p <- length(pips)

  x.diff <- get_hpm(object) - truth
  no.false.inclusions <- sum(x.diff > 0)
  no.false.exclusions <- sum(x.diff < 0)
  no.correct.choices <- sum(x.diff == 0)

  res <- list(false.inc = no.false.inclusions, false.exc = no.false.exclusions,
              false = p - no.correct.choices, brier = brier, p = p)
  class(res) <- "ipriorBVS_predict"
  res
}

#' @export
print.ipriorBVS_predict <- function(x) {
  cat("Brier score =", x$brier, "\n")
  cat("False choices =", x$false, "out of", x$p)
  cat(paste0(" (", dec_plac(x$false / x$p * 100) , "%)\n"))
  cat("\nOf these,\n")
  cat("False inclusions =", x$false.inc, "\n")
  cat("False exclusions =", x$false.exc, "\n")
}

#' @export
get_predict <- function(object, newdata, y.test = NULL) {
  beta <- coef(object)$tab[, "Mean"]
  X <- newdata
  y.hat <- beta[1] + X %*% beta[-1]
  if (!is.null(y.test)) test.error <- mean((y.test - y.hat) ^ 2)
  else test.error <- NULL

  list(y.hat = as.numeric(y.hat), mse = test.error)
}

#' @export
get_R2 <- function(object) {
  model.var <- 1 / summary(object)["psi", "Mean"]
  y.var <- attr(object$y, "scaled:scale")
  if (is.null(y.var)) y.var <- 1
  else y.var <- y.var ^ 2
  res <- 1 - model.var / y.var
  res
}