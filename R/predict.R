#' @export
predict.ipriorBVS <- function(object, truth = 1, ...) {
  pips <- get_pips(object)
  brier <- mean((pips - truth) ^ 2)
  brier
}