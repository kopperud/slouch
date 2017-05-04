bias_correction <- function(beta1, beta1.var, X, V, which.fixed.cov, which.random.cov, tree, pars, control, seed){
  list2env(tree, envir = environment())
  list2env(pars, envir = environment())
  list2env(control, envir = environment())
  list2env(seed, envir = environment())
  
  ## ###################################################################### ##
  ##                                                                        ##
  ##                            Bias correction                             ##
  ##            See equation (A.5) in Hansen & Bartoszek (2012)             ##
  ##                                                                        ##
  ## ###################################################################### ##
  
  n <- nrow(X)
  
  ## Vu and Vd are square block-matrices of size (n*ncol(X)).
  Vu <- diag(c(rep(0, n*length(beta1[-c(which.fixed.cov, which.random.cov),])), c(me.fixed.cov, me.random.cov)))
  Vd <- diag(c(rep(0, n*length(beta1[-c(which.fixed.cov, which.random.cov),])), c(sapply(Vd, function(e) diag(e)))))
  
  a_intercept <- X[, -c(which.fixed.cov, which.random.cov), drop = FALSE]
  
  if(!is.null(fixed.cov)){
    a_fixed <- matrix(apply(X[,which.fixed.cov, drop = FALSE], 
                            2, 
                            function(e) mean(e)), 
                      byrow = TRUE, nrow = nrow(X), 
                      ncol = length(which.fixed.cov))
  }else{
    a_fixed <- NULL
  }
  
  if(!is.null(random.cov)){
    a_random <- matrix(theta.X, ncol=length(theta.X), nrow = nrow(X), byrow = TRUE)
  }else{
    a_random <- NULL
  }
  
  a_matrix <- cbind(a_intercept, a_fixed, a_random)
  correction <- matrix(Vu%*%pseudoinverse(Vd+Vu)%*%c(X - a_matrix),  ncol=ncol(X), nrow=nrow(X), byrow=F)
  bias_corr <- pseudoinverse(t(X)%*%solve(V, X))%*%t(X)%*%solve(V, correction)
  m<-length(beta1)
  K <- solve(diag(1,m,m)-bias_corr)
  
  coefficients_bias_corr <- K%*%beta1
  #coefficients_var_bias_corr <- solve(K)%*%beta1.var%*%t(pseudoinverse(K))
  coefficients_var_bias_corr <- solve(K)%*%beta1.var%*%t(solve(K))
  
  res <- list(coefficients_bias_corr = matrix(cbind(coefficients_bias_corr, 
                                                    sqrt(diag(coefficients_var_bias_corr))),
                                              nrow = ncol(X),
                                              dimnames = list(colnames(X), c("Estimates", "Std. error"))),
              residuals_bias_corr = Y - (X %*% coefficients_bias_corr),
              K = K)
  return(res)
}