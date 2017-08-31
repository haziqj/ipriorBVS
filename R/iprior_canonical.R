iprior_canonical <- function(y, X) {
  p <- ncol(X)
  n <- nrow(X)
  psi <- 1
  lambda <- abs(rnorm(p))
  theta <- log(c(lambda, psi))
  y <- scale(y, scale = FALSE)
  X <- scale(X, scale = FALSE)
  XTX <- crossprod(X, X)

  environment(canonical_likelihood) <- environment()
  res <- optim(theta, canonical_likelihood, method = "BFGS", y = y, XTX = XTX,
               likelihood = TRUE, control = list(fnscale = -2))

  theta <- res$par
  res2 <- canonical_likelihood(theta, y, XTX, likelihood = FALSE)
  lambda <- res2$lambda
  psi <- res2$psi
  Vy.inv <- res2$Vy.inv
  LXt <- t(X) * lambda
  fitted.values <- as.numeric(psi * t(LXt) %*% (XTX %*% (LXt %*% (Vy.inv %*% y)))) +
    attr(y, "scaled:center")
  list(lambda = lambda, psi = psi, fitted.values = fitted.values, loglik = res$val)
}

canonical_likelihood <- function(theta, y, XTX, likelihood = TRUE) {
  p <- ncol(XTX)
  n <- length(y)
  lambda <- exp(theta[-(p + 1)])
  psi <- exp(theta[p + 1])
  psi.LXTXL <- (psi ^ 2) * (XTX * lambda) * rep(lambda, each = p)
  psi.LXTXL.inv <- chol2inv(chol(psi.LXTXL))
  A <- chol2inv(chol(psi.LXTXL.inv + XTX)) %*% t(X)
  # A <- solve(psi.LXTXL.inv + XTX, t(X))
  B <- X %*% A
  diag(B) <- diag(B) - 1
  Vy.inv <- -psi * B
  l <- (-n / 2) * log(2 * pi) + (1 / 2) * determinant(Vy.inv)$mod - 0.5 * crossprod(y, (Vy.inv %*% y))

  if (isTRUE(likelihood)) return(as.numeric(l))
  else return(list(
    lambda = lambda, psi = psi, Vy.inv = Vy.inv
  ))
}

canonical_old_school <- function(theta, y, X, likelihood = TRUE) {
  p <- ncol(X)
  n <- length(y)
  lambda <- exp(theta[-(p + 1)])
  psi <- exp(theta[p + 1])
  Xl <- split(X, rep(1:ncol(X), each = nrow(X)))
  Hl <- lapply(Xl, iprior::fnH2)
  H.lam <- Reduce("+", mapply(Hl, lambda, FUN = "*", SIMPLIFY = FALSE))
  class(H.lam) <- NULL

  tmp <- iprior::eigenCpp(H.lam)
  u <- tmp$val
  V <- tmp$vec


  Vy.inv.y <- as.numeric(crossprod(y, V) %*% (t(V) / (psi * u ^ 2 + 1 / psi)))
  l <- (-n / 2) * log(2 * pi) - (1 / 2) * sum(log(psi * u ^ 2 + 1 / psi)) - 0.5 * crossprod(y, Vy.inv.y)

  if (isTRUE(likelihood)) return(as.numeric(l))
  else return(list(
    lambda = lambda, psi = psi, Vy.inv = V %*% (t(V) / (psi * u ^ 2 + 1 / psi))
  ))

}