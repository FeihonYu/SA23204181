# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' @title A computer of Ranked Probability Score
#' @param pred_probs the prediction probability matrix
#' @param obs_outcomes the match outcome matrix
#' @export
RPS <- function(pred_probs, obs_outcomes) {
    .Call('_SA23204181_RPS', PACKAGE = 'SA23204181', pred_probs, obs_outcomes)
}

#' @title A computer of prediction accuracy
#' @param pred_probs the prediction probability matrix
#' @param obs_result the winner index vector
#' @export
accuracy <- function(pred_probs, obs_result) {
    .Call('_SA23204181_accuracy', PACKAGE = 'SA23204181', pred_probs, obs_result)
}

#' @title A computer of prediction accuracy
#' @param pred_probs the prediction probability matrix
#' @export
log_sum_max <- function(pred_probs) {
    .Call('_SA23204181_log_sum_max', PACKAGE = 'SA23204181', pred_probs)
}

#' @title A computer of prediction accuracy
#' @param pred_probs the prediction probability matrix
#' @export
sum_max <- function(pred_probs) {
    .Call('_SA23204181_sum_max', PACKAGE = 'SA23204181', pred_probs)
}

#' @title A computer of prediction accuracy
#' @param pred_probs the prediction probability matrix
#' @param odds the profit odds given by gamble company
#' @export
ideal_profit_sum <- function(pred_probs, odds) {
    .Call('_SA23204181_ideal_profit_sum', PACKAGE = 'SA23204181', pred_probs, odds)
}

#' @title A computer of prediction accuracy
#' @param pred_probs the prediction probability matrix
#' @param odds the profit odds given by gamble company
#' @param obs_outcomes the match outcome matrix
#' @export
actual_profit_sum <- function(pred_probs, odds, obs_outcomes) {
    .Call('_SA23204181_actual_profit_sum', PACKAGE = 'SA23204181', pred_probs, odds, obs_outcomes)
}

#' @title A computer of prediction accuracy
#' @param newdata the validation set of match outcomes 
#' @export
ComputeOut <- function(newdata) {
    .Call('_SA23204181_ComputeOut', PACKAGE = 'SA23204181', newdata)
}

#' @title A Gibbs sampler using Rcpp
#' @description A Gibbs sampler using Rcpp
#' @param N the number of samples
#' @param thin the number of between-sample random numbers
#' @param n,a,b the shape parameters concerned with binomial and beta distribution
#' @return a random sample of size \code{n}
#' @examples
#' \dontrun{
#' gibbs_cpp(100, 10 , 10,5,5)
#' }
#' @export
gibbs_cpp <- function(N, thin, n, a, b) {
    .Call('_SA23204181_gibbs_cpp', PACKAGE = 'SA23204181', N, thin, n, a, b)
}

