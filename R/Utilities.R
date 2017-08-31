dec_plac <- function(x, k = 2) format(round(x, k), nsmall = k)

decPlac <- dec_plac

is.ipriorBVS <- function(x) inherits(x, "ipriorBVS")

is.ipriorBVS_data <- function(x) inherits(x, "ipriorBVS_data")

as.data.frame.ipriorBVS_data <- function(x) x$data