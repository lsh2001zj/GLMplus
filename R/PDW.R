#' @title Weighted Poisson Deviance
#'
#' @param pred predictions
#' @param obs observations
#' @param ex weights
#'
#' @return t
#' @export
#'
#'
PDW <- function(pred, obs, ex = rep(1, length(obs))) {
  200 * sum( ex * ( pred - obs  + log((obs / pred) ^ (obs )))) / sum(ex)
}
