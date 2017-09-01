index_correct_choices <- function(x){
  paste(x, collapse = "")
}

# combine_mcmc <- function(x) {
#   res <- suppressWarnings(runjags::combine.mcmc(x))
# }

#' @export
get_gam_dat <- function(x) {
  if (is.ipriorBVS(x)) x <- x$mcmc
  mod.fit.mcmc <- suppressWarnings(suppressWarnings(coda::as.mcmc(x)))
  ind <- grep("gamma", colnames(mod.fit.mcmc))
  mod.fit.mcmc[, ind]
}

#' @export
index_models <- function(x) {
  if (is.ipriorBVS(x)) x <- x$mcmc
  apply(get_gam_dat(x), 1, index_correct_choices)
}

#' @export
get_dev_dat <- function(x) {
  if (is.ipriorBVS(x)) x <- x$mcmc
  mod.fit.mcmc <- suppressWarnings(coda::as.mcmc(x))
  ind <- grep("deviance", colnames(mod.fit.mcmc))
  mod.fit.mcmc[, ind]
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
get_model_names <- function(x) {
  if (is.ipriorBVS(x)) x <- x$mcmc
  models <- index_models(x)
  names(sort(table(models), decreasing = TRUE))
}

#' @export
get_pips <- function(x) {
  res <- apply(get_gam_dat(x), 2, mean)
  if (is.ipriorBVS(x)) names(res) <- x$xname
  res
}

#' @export
get_deviances <- function(x) {
  if (is.ipriorBVS(x)) x <- x$mcmc
  deviances <- get_dev_dat(x)
  models.list <- index_models(x)
  res <- unique.models <- get_model_names(x)
  for (i in seq_along(unique.models)) {
    res[i] <- mean(deviances[models.list == unique.models[i]])
  }
  res <- as.numeric(res)
  names(res) <- unique.models
  res
}

#' @export
get_pmps <- function(x) {
  if (is.ipriorBVS(x)) x <- x$mcmc
  models.list <- index_models(x)
  sort(table(models.list)  / length(models.list), decreasing = TRUE)
}

#' @export
get_hpm <- function(x) {
  res <- get_pmps(x)[1]
  res <- as.numeric(unlist(strsplit(names(res), split = "")))
  if (is.ipriorBVS(x)) names(res) <- x$xname
  res
}

#' @export
get_mpm <- function(x) {
  res <- tmp <- get_pips(x)
  res[] <- 0
  res[tmp >= 0.5] <- 1
  res
}

#' @export
tabulate_models <- function(x) {
  if (is.ipriorBVS(x)) x <- x$mcmc
  tab <- get_pmps(x)
  res <- strsplit(names(tab), split = "")
  res <- as.data.frame(lapply(res, as.numeric))
  res <- rbind(res, tab)
  res <- rbind(res, tab / tab[1])  # Bayes Factors
  res <- rbind(res, get_deviances(x))  # Deviances
  res <- cbind(c(get_pips(x), NA, NA, NA), res)  # Post. incl. probs.
  colnames(res) <- c("PIP", seq_len(ncol(res) - 1))
  prettify_table(res)
  # res
}

prettify_table <- function(x) {
  res <- x
  p <- nrow(res) - 3
  res[, 1] <- c(dec_plac(x[1:p, 1], k = 3), "", "", "")
  res[-nrow(res), -1] <- apply(x[-nrow(res), -1, drop = FALSE], 2, dec_plac, k = 3)
  res[nrow(res), -1] <- dec_plac(x[nrow(res), -1, drop = FALSE], k = 2)
  res[1:p, -1] <- apply(res[1:p, -1, drop = FALSE], 2, change_10)
  rownames(res)[-(1:p)] <- c("PMP", "BF", "Deviance")
  as.data.frame(res)
}

change_10 <- function(x) {
  res <- rep("x", length(x))
  res[grep("0.000", x)] <- ""
  res
}
