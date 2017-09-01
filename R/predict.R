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