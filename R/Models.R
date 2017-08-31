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
    gamma[j] ~ dbern(0.5)
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

#data# y, X, XTX.inv, n, p, pi
#inits# alpha, beta, gamma, psi, lambda
#monitor# gamma, beta, alpha, psi, deviance
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
  for (j in 1:p) { gamma[j] ~ dbern(0.5) }
  for (j in 1:p) { mu0[j] <- 0 }
  for (j in 1:p) { lambda[j] ~ dunif(0, 100) }
  for (j in 1:p) {
    for (k in 1:p) {
      lambda.inv[j,k] <- equals(j,k) / lambda[k]
    }
  }
  B <- lambda.inv[1:p, 1:p] %*% XTX.inv %*% lambda.inv[1:p, 1:p] / psi
  beta[1:p] ~ dmnorm(mu0, B)
}

#data# y, X, XTX.inv, n, p
#inits# alpha, beta, gamma, psi, lambda
#monitor# gamma, beta, alpha, psi, deviance
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
  for (j in 1:p) { gamma[j] ~ dbern(0.5) }
  for (j in 1:p) { mu0[j] <- 0 }
  beta[1:p] ~ dmnorm(mu0, B / psi)
}

#data# y, X, B, n, p
#inits# alpha, beta, gamma, psi
#monitor# gamma, beta, alpha, psi, deviance
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
    gamma[j] ~ dbern(0.5)
    beta[j] ~ dnorm(0, 0.001)
  }
  alpha ~ dnorm(0, 0.001)

  # Deviance
  for (i in 1:n) {
    d[i] <- log(2 * pi) - log(psi) + psi * pow(y[i] - mu[i], 2)
  }
  deviance <- sum(d)
}

#data# y, X, n, p, pi
#inits# alpha, beta, gamma, psi
#monitor# gamma, beta, alpha, psi, deviance
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
    gamma[j] ~ dbern(0.5)
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

#data# y, X, n, p, pi
#inits# alpha, beta, gamma, psi
#monitor# gamma, beta, alpha, psi, deviance
"

