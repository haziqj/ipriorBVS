#' @export
ipriorBVS <- function(...) {
  UseMethod("ipriorBVS")
}

#' @export
ipriorBVS.default <- function(y = Y, x = X, mod = mod.ipriorBVS,
                              control = list(), ...) {
  X <- as.matrix(x)
  Y <- as.vector(y)
  XTX <- crossprod(X)
  XTX.inv <- solve(XTX)
  p <- ncol(X)
  n <- nrow(X)

  # Default control for jags ---------------------------------------------------
  con <- list(n.chains = 4, n.iter = 5000, n.burnin = 2500, n.thin = 1)
  con_names <- names(con)
  con[(control_names <- names(control))] <- control
  if (length(noNms <- control_names[!control_names %in% con_names])) {
    warning("Unknown names in control options: ", paste(noNms, collapse = ", "),
            call. = FALSE)
  }
  control <- con

  # Linear model ---------------------------------------------------------------
  mod.lm <- lm(Y ~ 1 + X)
  beta.ols <- coef(summary(mod.lm))[, 1]
  sigma.ols <- summary(mod.lm)$sigma
  sigma.beta <- summary(mod.lm)$coefficients[, 2]

  # Bayesian Variable Selection ------------------------------------------------
  mod.data <- list(Y = Y, X = X, XTX.inv = XTX.inv, n = n, p = p)
  mod.params <- c("alpha", "Beta", "Gamma")
  mod.inits <- function() {
    list(tau = 1 / sigma.ols ^ 2, Beta = beta.ols[-1], Gamma = rep(1,p),
         lambda = 1)
  }
  mod.fit <- R2jags::jags(
    data  = mod.data,
    inits = mod.inits,
    parameters.to.save = mod.params,
    n.chains   = control$n.chains,
    n.iter     = control$n.iter,
    n.burnin   = control$n.burnin,
    model.file = mod,
    n.thin     = control$n.thin,
    DIC        = FALSE
  )

  # Results --------------------------------------------------------------------
  mod.fit.mcmc <- coda::as.mcmc(mod.fit)
  gam.ind <- grep("Gam", colnames(runjags::combine.mcmc(mod.fit.mcmc)))
  gam.dat <- runjags::combine.mcmc(mod.fit.mcmc)[, gam.ind]
  # res <- fn.resy(gam.dat)
  res <- list(mod.fit = mod.fit, gam.dat = gam.dat, xnames = paste0("X", 1:p))

  # Output ---------------------------------------------------------------------
  class(res) <- "ipriorBVS"
  res
}

ipriorBVS.formula <- function(formula, data = parent.frame(),
                              mod = mod.ipriorBVS, control = list(), ...) {
  mf <- model.frame(formula = formula, data = data)
  tt <- terms(mf)
  Terms <- delete.response(tt)
  X <- model.frame(Terms, mf)
  Y <- model.response(mf)
  res <- ipriorBVS.default(y = Y, x = X, mod = mod, control = control)
  res
}

#' @export
print.ipriorBVS <- function(x) {
  gam.summary <- fn.resy(x$gam.dat, x$xnames)
  print(head(gam.summary, 5))
}

#' @export
summary.ipriorBVS <- function(x) {
  print(x$mod.fit)
}

## I-prior BVS -----------------------------------------------------------------
#' @export
mod.ipriorBVS <- function() {
  for (j in 1:p) { gb[j] <- Gamma[j] * Beta[j] }

  for (i in 1:n) {
    Y[i] ~ dnorm(mu[i], tau)
    mu[i] <- alpha + inprod(X[i,1:p], gb[1:p])
  }

  tau ~ dgamma(0.01, 0.01)
  for (j in 1:p) { Gamma[j] ~ dbern(0.5) }
  Beta[1:p] ~ dmnorm(mu0, XTX.inv / (tau * lambda ^ 2))
  alpha ~ dnorm(0, 0.01)
  for (j in 1:p) { mu0[j] <- 0 }
  lambda ~ dunif(0, 100)
}

## Functions to analyse results ------------------------------------------------
fn4a <- function(x, varname){ 	#picks out correct choices.
  y <- as.numeric(x)
  if (sum(y) == 0) tmp <- "Intercept only"
  else tmp <- paste((varname)[y == 1], sep = "", collapse = " ")
  tmp
}

# fn4b <- function(x){		#counts false choices
# y <- as.numeric(x)
# this <- as.numeric(beta.true > 0)
# tmp <- as.character(sum(y != this))
# if(tmp == "0") tmp <- "None"
# tmp
# }

# fn4c <- function(x){		#counts false choices
# y <- as.numeric(x)
# this <- as.numeric(beta.true > 0)
# tmp <- as.character(sum(y != this))
# tmp
# }

#' @export
fn.resy <- function(x, xnames){
  model <- apply(x, 1, fn4a, varname = xnames)
  tab <- sort(table(model)  / length(model), decreasing = TRUE)
  ind <- (1:length(tab))[row.names(tab) == "None"]
  if (sum(row.names(tab) == "None") == 0) ind <- NA
  tab <- cbind(tab, tab[ind] / tab)
  colnames(tab) <- c("prop", "odds")
  # tab[order(rownames(tab)), ]
  tab
}

