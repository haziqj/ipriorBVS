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

  dat <- list(data = data.frame(y = y, X = X), truth = beta.true, snr = snr,
              sigma = sd)
  class(dat) <- "ipriorBVS_data"
  dat
}

#' @export
print.ipriorBVS_data <- function(x) {
  cat("n    = ", nrow(x$data), "\n")
  cat("p    = ", length(x$truth), "\n")
  if (!is.null(x$snr)) cat("SNR  = ", x$snr, "\n")
  beta <- matrix(x$truth)
  Sigma <- var(x$data[, -1])
  sigma <- x$sigma
  cat("SNR2 = ", as.numeric((t(beta) %*% Sigma %*% beta) / (sigma ^ 2)), "\n")
}

#' @export
gen_benchmark <- function(n = 40, sd = 1, seed = NULL) {
  p <- 8
  beta.true <- c(3, 1.5, 0, 0, 2, 0, 0, 0)
  if (!is.null(seed)) set.seed(seed)
  Sigma <- matrix(1, ncol = p, nrow = p)
  for (i in seq_len(p)) {
    for (j in seq_len(p)) {
      Sigma[i, j] <- 0.5 ^ abs(i - j)
    }
  }

  X <- mvtnorm::rmvnorm(n, mean = rep(0, p), sigma = Sigma)
  y <- X %*% matrix(beta.true) + rnorm(n, mean = 0, sd = sd)

  dat <- list(data = data.frame(y = y, X = X), truth = beta.true, sigma = sd)
  class(dat) <- "ipriorBVS_data"
  dat
}