# include <Rcpp.h>
using namespace Rcpp;
//' @title A computer of Ranked Probability Score
//' @param pred_probs the prediction probability matrix
//' @param obs_outcomes the match outcome matrix
//' @export
// [[Rcpp::export]]
double RPS(NumericMatrix pred_probs, IntegerMatrix obs_outcomes) {
  int num_matches = pred_probs.nrow();
  int num_categories = pred_probs.ncol();
  
  NumericVector rps_scores(num_matches);
  
  for (int i = 0; i < num_matches; ++i) {
    NumericVector pred_prob = pred_probs(i, _);
    IntegerVector obs_outcome = obs_outcomes(i, _);
    
    double rps = 0;
    
    NumericVector cum_probs = cumsum(pred_prob);
    NumericVector cum_outcomes = cumsum(as<NumericVector>(obs_outcome));
    
    for (int j = 0; j < num_categories; ++j) {
      rps += pow(cum_probs[j] - cum_outcomes[j], 2);
    }
    
    rps_scores[i] = rps / (num_categories - 1);
  }
  
  double mean_rps = mean(rps_scores);
  return mean_rps;
}
//' @title A computer of prediction accuracy
//' @param pred_probs the prediction probability matrix
//' @param obs_result the winner index vector
//' @export
// [[Rcpp::export]]
double accuracy(NumericMatrix pred_probs, IntegerVector obs_result) {
  int num_matches = pred_probs.nrow();
  int correct_predictions = 0;
  
  for (int i = 0; i < num_matches; ++i) {
    NumericVector pred_prob = pred_probs(i, _);
    int pred_index = which_max(pred_prob);

    if ((pred_index+1) == obs_result[i]) {
      correct_predictions += 1;
    }
  }
  
  double accuracy = (double) correct_predictions / num_matches;
  return accuracy;
}
//' @title A computer of prediction accuracy
//' @param pred_probs the prediction probability matrix
//' @export
// [[Rcpp::export]]
double log_sum_max(NumericMatrix pred_probs) {
  int num_prec = pred_probs.nrow();
  NumericVector max_probs(num_prec);
  
  for (int i = 0; i < num_prec; ++i) {
    max_probs[i] = max(pred_probs(i, _));
  }
  
  double log_sum = sum(log(max_probs));
  return log_sum;
}
//' @title A computer of prediction accuracy
//' @param pred_probs the prediction probability matrix
//' @export
// [[Rcpp::export]]
double sum_max(NumericMatrix pred_probs) {
  int num_prec = pred_probs.nrow();
  NumericVector max_probs(num_prec);
  
  for (int i = 0; i < num_prec; ++i) {
    max_probs[i] = max(pred_probs(i, _));
  }
  
  double sum_like = sum(max_probs);
  return sum_like;
}
//' @title A computer of prediction accuracy
//' @param pred_probs the prediction probability matrix
//' @param odds the profit odds given by gamble company
//' @export
// [[Rcpp::export]]
double ideal_profit_sum(NumericMatrix pred_probs, NumericVector odds) {
  int num_prec = pred_probs.nrow();
  double expec_profit = 0;
  IntegerVector index(num_prec);
  
  for (int i = 0; i < num_prec; ++i) {
    index[i] = which_max(pred_probs(i, _));
    expec_profit += pred_probs(i, index[i]) * odds[i] - 1;
  }
  
  return expec_profit / num_prec;
}
//' @title A computer of prediction accuracy
//' @param pred_probs the prediction probability matrix
//' @param odds the profit odds given by gamble company
//' @param obs_outcomes the match outcome matrix
//' @export
// [[Rcpp::export]]
double actual_profit_sum(NumericMatrix pred_probs, NumericVector odds, IntegerMatrix obs_outcomes) {
  int num_prec = pred_probs.nrow();
  double expec_profit = 0;
  IntegerVector index(num_prec);
  
  for (int i = 0; i < num_prec; ++i) {
    index[i] = which_max(pred_probs(i, _));
    
    if (obs_outcomes(i, index[i]) == 1) 
      expec_profit += (odds[i] - 1);
    else 
      expec_profit -= 1;
  }
  
  return expec_profit / num_prec;
}
//' @title A computer of prediction accuracy
//' @param newdata the validation set of match outcomes 
//' @export
// [[Rcpp::export]]
NumericMatrix ComputeOut(DataFrame newdata) {
  CharacterVector outcomes = newdata["FTR"];
  IntegerVector levels = match(outcomes, CharacterVector::create("A", "D", "H"));
  
  int numRows = newdata.nrows();
  NumericMatrix result(numRows, 3);
  
  for (int i = 0; i < numRows; ++i) {
    result(i, levels[i] - 1) = 1;
  }
  
  return result;
}  
//' @title A Gibbs sampler using Rcpp
//' @description A Gibbs sampler using Rcpp
//' @param N the number of samples
//' @param thin the number of between-sample random numbers
//' @param n,a,b the shape parameters concerned with binomial and beta distribution
//' @return a random sample of size \code{n}
//' @examples
//' \dontrun{
//' gibbs_cpp(100, 10 , 10,5,5)
//' }
//' @export
// [[Rcpp::export]]
 NumericMatrix gibbs_cpp(int N, int thin ,int n, double a ,double b) {
    NumericMatrix mat(N, 2);
    double x = 0, y = 0;
    for(int i = 0; i < N; i++) {
      for(int j = 0; j < thin; j++) {
        x = rbinom(1, n, y)[0];
        y = rbeta(1,x+a,n-x+b)[0];
      }
      mat(i, 0) = x;
      mat(i, 1) = y;
    }
    return(mat);
  }