bvs_iprior_sing <-
"model{
  for (j in 1:p) {
    gb[j] <- gamma[j] * beta[j]
  }

  for (i in 1:n) {
    y[i] ~ dnorm(mu[i], psi)
    mu[i] <- alpha + inprod(X[i, 1:p], gb[1:p])
  }

  # Priors
  psi ~ dgamma(0.001, 0.001)
  for (j in 1:p) {
    gamma[j] ~ dbern(gamma.prob[j])
  }
  beta[1:p] ~ dmnorm(mu0, XTX.inv / (psi * lambda ^ 2))
  alpha ~ dnorm(0, 0.001)
  for (j in 1:p) {
    mu0[j] <- 0
  }
  lambda ~ dgamma(0.001, 0.001)

  # Deviance
  for (i in 1:n) {
    d[i] <- log(2 * pi) - log(psi) + psi * pow(y[i] - mu[i], 2)
  }
  deviance <- sum(d)
}

#data# y, X, XTX.inv, n, p, pi, gamma.prob
#inits# alpha, beta, gamma, psi, lambda
#monitor# gamma, alpha, gb, psi, deviance
"

bvs_iprior_mult <- "model{
  for (j in 1:p) {
    gb[j] <- gamma[j] * beta[j]
  }

  for (i in 1:n) {
    y[i] ~ dnorm(mu[i], psi)
    mu[i] <- alpha + inprod(X[i, 1:p], gb[1:p])
  }

  # Priors
  psi ~ dgamma(0.001, 0.001)
  alpha ~ dnorm(0, 0.001)
  for (j in 1:p) {
    gamma[j] ~ dbern(gamma.prob[j])
    mu0[j] <- 0
    lambda[j] ~ dgamma(0.1, 0.1)
  }
  for (j in 1:p) {
    for (k in 1:p) {
      lambda.inv[j,k] <- equals(j,k) / lambda[k]
    }
  }
  B <- lambda.inv[1:p, 1:p] %*% XTX.inv %*% lambda.inv[1:p, 1:p] / psi
  beta[1:p] ~ dmnorm(mu0, B)

  # Deviance
  for (i in 1:n) {
    d[i] <- log(2 * pi) - log(psi) + psi * pow(y[i] - mu[i], 2)
  }
  deviance <- sum(d)
}

#data# y, X, XTX.inv, n, p, pi, gamma.prob
#inits# alpha, beta, gamma, psi, lambda
#monitor# gamma, alpha, gb, psi, deviance, lambda
"

bvs_iprior_mult_fixed <- "model{
  for (j in 1:p) {
    gb[j] <- gamma[j] * beta[j]
  }

  for (i in 1:n) {
    y[i] ~ dnorm(mu[i], psi)
    mu[i] <- alpha + inprod(X[i, 1:p], gb[1:p])
  }

  # Priors
  psi ~ dgamma(0.001, 0.001)
  alpha ~ dnorm(0, 0.001)
  for (j in 1:p) {
    gamma[j] ~ dbern(gamma.prob[j])
    mu0[j] <- 0
  }
  beta[1:p] ~ dmnorm(mu0, B / psi)

  # Deviance
  for (i in 1:n) {
    d[i] <- log(2 * pi) - log(psi) + psi * pow(y[i] - mu[i], 2)
  }
  deviance <- sum(d)
}

#data# y, X, B, n, p, pi, gamma.prob
#inits# alpha, beta, gamma, psi
#monitor# gamma, alpha, gb, psi, deviance
"

bvs_independent <- "model{
  for (j in 1:p) {
    gb[j] <- gamma[j] * beta[j]
  }

  for (i in 1:n) {
    y[i] ~ dnorm(mu[i], psi)
    mu[i] <- alpha + inprod(X[i, 1:p], gb[1:p])
  }

  # Priors
  psi ~ dgamma(0.001, 0.001)
  for (j in 1:p) {
    gamma[j] ~ dbern(gamma.prob[j])
    beta[j] ~ dnorm(0, 0.001)
  }
  alpha ~ dnorm(0, 0.001)

  # Deviance
  for (i in 1:n) {
    d[i] <- log(2 * pi) - log(psi) + psi * pow(y[i] - mu[i], 2)
  }
  deviance <- sum(d)
}

#data# y, X, n, p, pi, gamma.prob
#inits# alpha, beta, gamma, psi
#monitor# gamma, alpha, gb, psi, deviance
"

bvs_gprior <- "model{
  for (j in 1:p) {
    gb[j] <- gamma[j] * beta[j]
  }

  for (i in 1:n) {
    y[i] ~ dnorm(mu[i], psi)
    mu[i] <- alpha + inprod(X[i, 1:p], gb[1:p])
  }

  # Priors
  psi ~ dgamma(0.001, 0.001)
  for (j in 1:p) {
    gamma[j] ~ dbern(gamma.prob[j])
    mu0[j] <- 0
  }
  alpha ~ dnorm(0, 0.001)
  Beta[1:p] ~ dmnorm(mu0, XTX / n)

  # Deviance
  for (i in 1:n) {
    d[i] <- log(2 * pi) - log(psi) + psi * pow(y[i] - mu[i], 2)
  }
  deviance <- sum(d)
}

#data# y, X, n, p, pi, gamma.prob
#inits# alpha, beta, gamma, psi
#monitor# gamma, alpha, gb, psi, deviance
"


bvs_iprior2 <- "model{
  for (j in 1:p) { gl[j] <- gamma[j] * lambda[j] }

  # Mean vector
  for (i in 1:n) { mu[i] <- alpha }

  # Psi inverse matrix = diag(1 / psi)
  for (i in 1:n) {
    for (j in 1:n) {
      Psi.inv[i, j] <- equals(i, j) * psi
    }
  }

  # H.lam matrix
  for (j in 1:p) {
    Hlamj[1:n, 1:n, j] <- gl[j] * H[1:n, 1:n, j]
  }
  Hlam <- Psi.inv
  # for (i in 1:n) {
  #   for (j in 1:n) {
  #     for (k in 1:p) {
  #       Hlam[i,j] <- Hlamj[k,i,j]
  #     }
  #   }
  # }

  # Covariance matrix for marginal y
  Vy <- psi * Hlam %*% Hlam + Psi.inv

  y[1:n] ~ dmnorm(mu, Vy)

  # Priors
  alpha ~ dnorm(0, 0.001)
  psi ~ dgamma(0.001, 0.001)
  for (j in 1:p) { gamma[j] ~ dbern(0.5) }
  for (j in 1:p) { lambda[j] ~ dunif(0, 100) }

  # Deviance
  for (i in 1:n) {
    d[i] <- log(2 * pi) - log(psi) + psi * pow(y[i] - mu[i], 2)
  }
  deviance <- sum(d)
}

#data# y, H, n, p, pi
#inits# alpha, gamma, psi, lambda
#monitor# gamma, alpha, gb, psi, deviance
"
