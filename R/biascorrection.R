bias_correction <- function(beta1, beta1.var, Y, X, V, w_beta, c_regimes, seed, observations, control){
  ## ###################################################################### ##
  ##                                                                        ##
  ##                            Bias correction                             ##
  ##            See equation (A.5) in Hansen & Bartoszek (2012)             ##
  ##                                                                        ##
  ## ###################################################################### ##
  
  n <- nrow(X)
  
  ## Redefine Vu and Vd as square block-matrices of size (n*ncol(X)).
  # ncov <- length(beta1)# - ncol(observations$direct.cov) - ncol(observations$random.cov)
  # for (x in c("direct.cov", "random.cov")){
  #   if (!is.null(observations[[x]])){
  #     ncov <- ncov - ncol(observations[[x]])
  #   }
  # }
  
  if (control$interactions & !is.null(w_beta$which.direct.cov) & !is.null(w_beta$which.regimes)){
    Vu_direct <- seed$Vu[1:(length(w_beta$which.direct.cov) / length(w_beta$which.regimes))]
    Vu_direct <- do.call(c, 
                  lapply(head(seed$Vu, length(w_beta$which.direct.cov) / length(w_beta$which.regimes)), 
                         function(Vu_i) lapply(split(c_regimes, colnames(c_regimes)), 
                                               function(c_regime) c_regime*Vu_i)))
    Vu <- c(Vu_direct, 
            tail(seed$Vu, length(w_beta$which.random.cov)))
    
    Vd_direct <- do.call(c, 
                  lapply(head(seed$Vd, length(w_beta$which.direct.cov) / length(w_beta$which.regimes)), 
                         function(Vd_i) lapply(split(c_regimes, colnames(c_regimes)), 
                                               function(c_regime) c_regime*Vd_i)))
    Vd <- c(Vd_direct, 
            tail(seed$Vd, length(w_beta$which.random.cov)))
    
  }else{
    Vu <- seed$Vu
    Vd <- seed$Vd
  }
  
  zeros <- rep(0, n*(length(beta1) - length(beta1[c(w_beta$which.direct.cov, w_beta$which.random.cov)])))
  Vu <- diag(c(zeros, sapply(Vu, diag)))
  Vd <- diag(c(zeros, sapply(Vd, diag)))
  
  a_intercept <- X[, -c(w_beta$which.direct.cov, w_beta$which.random.cov), drop = FALSE]
  
  if(length(w_beta$which.direct.cov) > 0){
    a_direct <- matrix(apply(X[,w_beta$which.direct.cov, drop = FALSE], 
                            2, 
                            function(e) mean(e)), 
                      byrow = TRUE, nrow = nrow(X), 
                      ncol = length(w_beta$which.direct.cov))
  }else{
    a_direct <- NULL
  }
  
  if(length(w_beta$which.random.cov) > 0){
    a_random <- matrix(seed$brownian_mean, ncol=length(seed$brownian_mean), nrow = nrow(X), byrow = TRUE)
  }else{
    a_random <- NULL
  }
  
  #stop()
  a_matrix <- cbind(a_intercept, a_direct, a_random)
  correction <- matrix(Vu%*%pseudoinverse(Vd+Vu)%*%c(X - a_matrix),  ncol=ncol(X), nrow=nrow(X), byrow=F)
  bias_corr <- pseudoinverse(t(X)%*%solve(V, X))%*%t(X)%*%solve(V, correction)
  m<-length(beta1)
  K <- solve(diag(1,m,m)-bias_corr)
  colnames(K) <- rownames(K) <- colnames(X)

  
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