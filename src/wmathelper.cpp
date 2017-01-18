#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
NumericVector regimeWeightsHelper(double a, NumericVector x){
  int n = x.size();
  NumericVector res(n);
  double tip1;
  double ti;
  
  res[0] = 1 - exp(-a*(x[0] - x[1]));
  
  for (int i=1; i < n; i++){
    tip1 = exp(-a*(x[0] - x[i + 1]));
    
    ti = exp(-a*(x[0] - x[i]));
    res[i] = ti - tip1;
    ti = tip1;
  }
  res[n-1] = exp(-a*x[0]);
  
  return res;
}

// 
// res <- c(vapply(head(seq_along(nt), -1),
//                 function(i) exp(-a*(nt[1] - nt[i])) - exp(-a*(nt[1] - nt[i + 1])), FUN.VALUE = 0),
//                 exp(-a*nt[1]) ## Theta0\Ya


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R

# nt <- m4$lineages[[34]]$nodes_time
# regimeWeights(0.3, nt)
# print(res5)
# print(microbenchmark(regimeWeights(0.3, nt)))
*/
