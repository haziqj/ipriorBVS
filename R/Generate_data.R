#' @export
gen_bvs <- function(n = 150, p = 100, snr = 0.9, sd = 2, seed = NULL) {
  p.on <- floor(p * snr)
  p.off <- p - p.on
  beta.true <- c(rep(1, p.on), rep(0, p.off))
  if (!is.null(seed)) set.seed(seed)
  X <- matrix(rnorm(n * p), nrow = n)
  Z <- rnorm(n)
  X <- X + matrix(Z, nrow = n, ncol = p)
  y <- X %*% matrix(beta.true) + rnorm(n, mean = 0, sd = sd)

  dat <- list(data = data.frame(y = y, X = X), truth = beta.true, snr = snr)
  class(dat) <- "ipriorBVS_data"
  dat
}

#' @export
print.ipriorBVS_data <- function(x) {
  cat("n   = ", nrow(x$data), "\n")
  cat("p   = ", length(x$truth), "\n")
  cat("SNR = ", x$snr)
}

