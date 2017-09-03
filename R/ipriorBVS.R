.onLoad <- function(libname, pkgname) {
  runjags::runjags.options(
    silent.runjags = TRUE,
    silent.jags = FALSE,
    summary.warning = FALSE,
    rng.warning = FALSE
  )
}

#' @export
ipriorBVS <- function(...) {
  UseMethod("ipriorBVS")
}

#' @export
ipriorBVS.default <- function(y, X, model = "iprior_sing", two.stage = FALSE,
                              stand.x = TRUE, stand.y = TRUE,
                              n.chains = parallel::detectCores(),
                              n.samp = 10000, n.burnin = 4000, n.adapt = 1000,
                              n.thin = 1, n.par = n.chains, priors = NULL) {
  y <- as.numeric(scale(y, scale = stand.y, center = stand.y))
  X <- scale(X, scale = stand.x, center = stand.x)
  XTX <- crossprod(X)
  XTX.inv <- solve(XTX)
  p <- ncol(X)
  n <- nrow(X)
  xnames <- colnames(X)

  model <- match.arg(model, c(
    "iprior_sing",
    "iprior_mult",
    "iprior_mult_fixed",
    "flat_prior",
    "gprior",
    "iprior2"
  ))

  # # Default control for jags ---------------------------------------------------
  # con <- list(n.chains = 4, n.iter = 5000, n.burnin = 2500, n.thin = 1)
  # con_names <- names(con)
  # con[(control_names <- names(control))] <- control
  # if (length(noNms <- control_names[!control_names %in% con_names])) {
  #   warning("Unknown names in control options: ", paste(noNms, collapse = ", "),
  #           call. = FALSE)
  # }
  # control <- con

  # Linear model ---------------------------------------------------------------
  mod.lm <- lm(y ~ 1 + X)
  beta.ols <- coef(summary(mod.lm))[, 1]
  sigma.ols <- summary(mod.lm)$sigma
  sigma.beta <- summary(mod.lm)$coefficients[, 2]

  # Bayesian Variable Selection ------------------------------------------------
  pi <- pi
  # Initial values
  alpha <- beta.ols[1]
  beta <- beta.ols[-1]
  psi <- 1 / sigma.ols ^ 2
  gamma <- rep(1, p)
  gamma.prob <- rep(0.5, p)
  if (model == "iprior_sing") {
    bvs_model <- bvs_iprior_sing
    lambda <- 1
  }
  if (model == "iprior_mult") {
    bvs_model <- bvs_iprior_mult
    # mod <- iprior_canonical(y, X)
    # lambda.mle <- mod$lambda
    lambda <- rep(1, p)
  }
  if (model == "iprior_mult_fixed") {
    bvs_model <- bvs_iprior_mult_fixed
    mod <- iprior_canonical(y, X)
    lambda <- mod$lambda
    B <- diag(1 / lambda) %*% XTX.inv %*% diag(1 / lambda)
  }
  if (model == "flat_prior") {
    bvs_model <- bvs_independent
  }
  if (model == "gprior") {
    bvs_model <- bvs_independent
  }
  if (model == "iprior2") {
    bvs_model <- bvs_iprior2
    H <- array(NA, dim = c(n, n, p))
    for (j in 1:p) H[1:n, 1:n, j] <- iprior::fnH2(X[, j])
    lambda <- rep(1, p)
  }

  if (!is.null(priors)) {
    for (i in seq_along(priors)) {
      position.of.tilde <- regexpr("~", priors[i])[1]
      find.this <- substr(priors[i], 1, position.of.tilde)
      bvs_model <- gsub(find.this, paste(priors[i], "#"), bvs_model)
    }
  }

  # Run model
  mod.fit <- runjags::run.jags(bvs_model, n.chains = n.chains, burnin = n.burnin,
                               adapt = n.adapt, sample = n.samp / n.chains,
                               thin = n.thin, method = "parallel", n.sims = n.par)


  cat("\n")

  if (isTRUE(two.stage)) {
    cat("Running second stage variable selection\n")
    gamma.prob <- get_mpm(mod.fit) * 0.5
    mod.fit <- runjags::run.jags(bvs_model, n.chains = n.chains, burnin = n.burnin,
                                 adapt = n.adapt, sample = n.samp / n.chains,
                                 thin = n.thin, method = "parallel", n.sims = n.par)
  }

  # Results --------------------------------------------------------------------
  res <- list(mcmc = mod.fit, xnames = xnames)

  # Output ---------------------------------------------------------------------
  class(res) <- "ipriorBVS"
  res
}

#' @export
ipriorBVS.formula <- function(formula, data = parent.frame(),
                              model = "iprior_sing", two.stage = FALSE,
                              stand.x = TRUE, stand.y = TRUE,
                              n.chains = parallel::detectCores(),
                              n.samp = 10000, n.burnin = 4000, n.adapt = 1000,
                              n.thin = 1, n.par = n.chains, priors = NULL, ...) {
  if (is.ipriorBVS_data(data)) data <- as.data.frame(data)
  mf <- model.frame(formula = formula, data = data)
  tt <- terms(mf)
  Terms <- delete.response(tt)
  X <- model.frame(Terms, mf)
  Y <- model.response(mf)
  colnames(X)
  res <- ipriorBVS.default(Y, X, model, two.stage, stand.x, stand.y, n.chains,
                           n.samp, n.burnin, n.adapt, n.thin, n.par, priors, ...)
  res
}

#' @export
print.ipriorBVS <- function(x, n = 5) {
  gam.summary <- tabulate_models(x$mcmc)
  rownames(gam.summary)[seq_along(x$xnames)] <- x$xnames
  no.of.cols <- min(n + 1, ncol(gam.summary))
  print(gam.summary[, seq_len(no.of.cols)])
}

#' @export
summary.ipriorBVS <- function(x) {
  summary(x$mcmc)
}




