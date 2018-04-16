bias_correction <- function(beta1, beta1.var, Y, X, V, which.direct.cov, which.random.cov, seed){
  ## ###################################################################### ##
  ##                                                                        ##
  ##                            Bias correction                             ##
  ##            See equation (A.5) in Hansen & Bartoszek (2012)             ##
  ##                                                                        ##
  ## ###################################################################### ##
  
  n <- nrow(X)
  
  ## Redefine Vu and Vd as square block-matrices of size (n*ncol(X)).
  zeros <- rep(0, n*length(beta1[-c(which.direct.cov, which.random.cov),]))
  Vu <- diag(c(zeros, sapply(seed$Vu, diag)))
  Vd <- diag(c(zeros, sapply(seed$Vd, diag)))
  
  a_intercept <- X[, -c(which.direct.cov, which.random.cov), drop = FALSE]
  
  if(length(which.direct.cov) > 0){
    a_fixed <- matrix(apply(X[,which.direct.cov, drop = FALSE], 
                            2, 
                            function(e) mean(e)), 
                      byrow = TRUE, nrow = nrow(X), 
                      ncol = length(which.direct.cov))
  }else{
    a_fixed <- NULL
  }
  
  if(length(which.random.cov) > 0){
    a_random <- matrix(seed$theta.X, ncol=length(seed$theta.X), nrow = nrow(X), byrow = TRUE)
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
  vcov_bias_corr <- solve(K)%*%beta1.var%*%t(solve(K))
  dimnames(vcov_bias_corr) <- list(colnames(X), colnames(X))
  
  res <- list(coefficients_bias_corr = matrix(cbind(coefficients_bias_corr, 
                                                    sqrt(diag(vcov_bias_corr))),
                                              nrow = ncol(X),
                                              dimnames = list(colnames(X), c("Estimates", "Std. error"))),
              residuals_bias_corr = Y - (X %*% coefficients_bias_corr),
              vcov_bias_corr = vcov_bias_corr,
              K = K)
  return(res)
}