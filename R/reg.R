## Calculate determinant of V
mk.log.det.V <- function(V, N){
  det.V<-det(V)
  if(det.V==0){
    print(paste("Warning: Determinant of V = 0"))
    #Minimum value of diagonal scaling factor
    inv.min.diag.V<-1/min(diag(V))
    V<-V*inv.min.diag.V
    #Rescale and log determinant
    log(det(V))+log(min(diag(V)))*N
  }
  else {
    log(det.V)
  }
}

## Design matrix
calc.X <- function(a, hl, treepar, modelpar, seed, is.opt.reg = TRUE){
  list2env(modelpar, envir = environment())
  list2env(treepar, envir = environment())
  list2env(seed, envir = environment())
  
  if(is.opt.reg == TRUE){
    rho <- (1-(1-exp(-a*T.term))/(a*T.term))
  }else{
    rho <- 1
  }
  if(!is.null(fixed.fact)){
    #m1 <- weight.matrix(a, topology, times, N, regime.specs, fixed.cov = NULL, intercept)
    m1 <- weight.matrix(a, ancestor, times, N, regime.specs, fixed.cov = NULL, intercept, weight.m.regimes = regimes1, ep = epochs1)
    #print(microbenchmark(weight.matrix(a, topology, times, N, regime.specs, fixed.cov = NULL, intercept, weight.m.regimes = regimes1, ep = epochs1)))
    if(!is.null(fixed.cov) | !is.null(random.cov)){
      m2 <- matrix(cbind(fixed.pred, rho*pred), nrow = N, dimnames = list(NULL, c(names.fixed.cov, names.random.cov)))
    }else{
      m2 <- NULL
    }
    return(cbind(m1,m2))
  }else{
    if(!is.null(modelpar$intercept)){
      K <- 1
      K.name <- "Intercept"
    }else{
      K <- cbind(exp(-a*T.term),
                 1-exp(-a*T.term),
                 if (!is.null(modelpar$random.cov) | !is.null(modelpar$fixed.cov)) 1-exp(-a*T.term)-(1-(1-exp(-a*T.term))/(a*T.term)) else NULL)
      K.name <- c("Ya",
                  "b0",
                  if (!is.null(modelpar$random.cov) | !is.null(modelpar$fixed.cov)) "X_anc" else NULL)
    }
    matrix(cbind(K,
                 fixed.pred,
                 rho*pred), 
           nrow=N,
           dimnames=list(NULL, c(K.name, names.fixed.cov, names.random.cov)))
  }
}

# Part "two" of variance-covariance matrix, hansen et al. 2008.
# It looks a little cryptic since it was vectorized from a nested for-loop
# Still is a bottleneck in performance. 01 sept 2016. Consider rewrite in Rcpp ?
calc.cm2 <- function(a, T.term, N, tia, tja, ta){
  ti <- matrix(rep(T.term, N), nrow=N)
  tj <- t(ti)
  term0 <- (1-exp(-a*ti))
  
  num.prob <- ifelse(ta == 0, 1, (1-exp(-a*ta))/(a*ta))
  return((term0/(a*ti))*(t(term0)/(a*tj)) -
           (exp(-a*tia)*t(term0)/(a*tj) + exp(-a*tja)*term0/(a*ti))*num.prob)
  # return((term0/(a*ti))*(t(term0)/(a*tj)) -  ##### INTENTIONALLY WRONG
  #          (exp(-a*tia)*(1 - exp(-a*ti))/(a*tj) + exp(-a*tja)*term0/(a*ti))*num.prob)
}

calc.V <- function(hl, vy, a, cm2, beta1, which.fixed.cov, which.random.cov, random.cov, T.term, fixed.cov, Vu_given_x, me.cov, n.pred, mecov.fixed.cov, n.fixed.pred, N, s.X, ta, tij, me.response){
  
  
  if (hl == 0 | is.null(random.cov)){
    rho2 <- 1
  }   else {
    rho2 <- (1-(1-exp(-a*T.term))/(a*T.term))^2
  }
  if(!is.null(fixed.cov) | !is.null(random.cov)){
    
    ## Update measurement error in X.
    ## Measurement error in predictor - take 2
    beta_continuous <- beta1[c(which.fixed.cov, which.random.cov)]
    obs_var_con2 <- list()
    for (i in 1:length(beta_continuous)){
      obs_var_con2[[i]] <- Vu_given_x[[i]]*(beta_continuous[i]^2)*rho2
    }
    obs_var_con2 <- Reduce('+', obs_var_con2)
    
    # ## BROKEN!? Needs test.
    # calculate covariances between response and the stochastic predictor, to be subtracted in the diagonal of V
    if(sum(me.cov) == 0){
      mcov <- 0
    }else{
      mcov <- diag(rowSums(matrix(data=as.numeric(me.cov)*t(kronecker(2*beta1[which.random.cov,],(1-(1-exp(-a*T.term))/(a*T.term)))), ncol=n.pred)))
    }
    
    if(sum(mecov.fixed.cov) == 0){
      mcov.fixed <- 0
    }else{
      mcov.fixed <- diag(rowSums(matrix(data=as.numeric(mecov.fixed.cov)*t(kronecker(2*beta1[which.fixed.cov,],(1-(1-exp(-a*T.term))/(a*T.term)))), ncol=n.fixed.pred)))
    }
  }else{
    obs_var_con2 <- 0
    mcov <- 0
    mcov.fixed <- 0
  }
  ## Piece together V
  if(hl == 0){
    cm0 <- diag(rep(vy, times=N))
  }else{
    if(!is.null(random.cov)){
      s1 <- sum(s.X%*%((beta1[which.random.cov,])^2))
      cm0 <- (s1/(2*a)+vy)*(1-exp(-2*a*ta))*exp(-a*tij) + (s1*ta*cm2)
    }else{
      cm0 <- vy*(1-exp(-2*a*ta))*exp(-a*tij)
    }
  }
  V <- cm0 + na.exclude(me.response) + obs_var_con2 - mcov - mcov.fixed
  return(V)
}


test.conv <- function(beta.i, beta1, con.count, convergence){
  # if(!is.null(modelpar$intercept)){
  #   test <- ifelse(abs(as.numeric(beta.i - beta1)) <= convergence, 0, 1)
  # }else if(is.null(modelpar$intercept) & (!is.null(modelpar$fixed.fact)) | (is.null(fixed.cov) & is.null(random.cov))){
  #   test <- ifelse(abs(as.numeric(beta.i - beta1))[-1] <= convergence, 0, 1)
  # }
  # else{
  #   test <- ifelse(abs(as.numeric(beta.i - beta1))[-(1:2)] <= convergence, 0, 1)
  #   ## Effectively removes beta[1:2] <= 0.001 from being criteria in convergence, when non-ultrametric with continuous covariates
  # }
  beta.i <- beta.i[!is.na(beta.i),]
  beta1 <- beta1[!is.na(beta1), ]
  test <- ifelse(abs(as.numeric(beta.i - beta1)) <= convergence, 0, 1)
  
  if(sum(test)==0) return (TRUE)
  if(con.count >= 50)
  {
    message("Warning, estimates did not converge after 50 iterations, last estimates printed out")
    return(TRUE)
  }
  return(FALSE)
}

reg <- function(hl_vy, modelpar, treepar, seed, gridsearch = TRUE){
  list2env(treepar, envir = environment())
  list2env(modelpar, envir = environment())
  list2env(seed, envir = environment())
  
  hl <- hl_vy[1]; vy <- hl_vy[2]
  
  if(hl == 0){
    a <- Inf ## When response is modeled only by fixed factors, of which one or more levels contain only internal regimes, Var(beta) becomes singular, throws error.
  }else{
    a <- log(2)/hl
    if (!is.null(random.cov)){
      cm2 <- calc.cm2(a, T.term, N, tia, tja, ta)
    }else cm2 <- NULL
  }
  X <- calc.X(a, hl, treepar, modelpar, seed, is.opt.reg = TRUE)
  
  ## ols.beta1 exists from OLS estimate
  #beta1 <- ols.beta1
  #beta1 <- solve(t(X)%*%X)%*%(t(X)%*%Y) ## Confirm: Is it ok to do qr-decomposition instead of usual beta estimator? Seems to avoid some bugs, e.g 
  ##  Error in solve.default(t(X) %*% X) : 
  ##  system is computationally singular: reciprocal condition number = 4.96821e-23 
  beta1 <- matrix(qr.coef(qr(X), Y), ncol=1)
  beta1.descriptor <- c(rep("Intercept", length(beta1) - length(names.fixed.cov) - length(names.random.cov)),
                        rep("Instantaneous cov", length(names.fixed.cov)),
                        rep("Random cov", length(names.random.cov)))
  which.fixed.cov <- which(beta1.descriptor == "Instantaneous cov")
  which.random.cov <- which(beta1.descriptor == "Random cov")
  
  con.count <- 0
  repeat{
    V <- calc.V(hl, vy, a, cm2, beta1, which.fixed.cov, which.random.cov, random.cov, T.term, fixed.cov, Vu_given_x, me.cov, n.pred, mecov.fixed.cov, n.fixed.pred, N, s.X, ta, tij, me.response)
    #print(microbenchmark(calc.V(hl, vy, a, cm2, beta1 = beta1)))
    
    
    #V.inverse <- pseudoinverse(V)
    
    if(FALSE){
      V.inverse <- solve(V)
      beta.i.var <- solve(t(X)%*%V.inverse%*%X)
      #beta.i.var <- pseudoinverse(t(X)%*%V.inverse%*%X)
      
      #### Ask Thomas about this one.
      if(Inf %in% beta.i.var) {print("Pseudoinverse of (XT * V * X) contained values = Inf, which were set to 10^300")}
      beta.i.var <-replace(beta.i.var, beta.i.var ==Inf, 10^300)
      if(-Inf %in% beta.i.var) {print("Pseudoinverse of (XT * V * X) contained values = -Inf, which were set to -10^300")}
      beta.i.var <-replace(beta.i.var, beta.i.var ==-Inf, -10^300)
      
      beta.i<-beta.i.var%*%(t(X)%*%V.inverse%*%Y)
      #if(all(V.inverse[!diag(nrow(V.inverse))] == 0)) warning("For hl = ", hl," and vy = ", vy," the inverse of V is strictly diagonal.")
    }else{
      # Linear transformation of both Y and model matrix, by the inverse of cholesky decomposition of V
      # Is this method ok to use? Seems to avoid some bugs (see below). R's lm.fit() does not throw
      # error when X is singular, but returns NA in the coefficients.
      # Problem with this transform is, does not work when matrix V is not positive semi-definite (Should this be possible? happens w dummy dataset on random.cov)
      #
      #  Error in solve.default(t(X) %*% V.inverse %*% X) : 
      #  system is computationally singular: reciprocal condition number = 2.28389e-22
      #
      # if(any(Re(eigen(V)$values) < 0)){
      #   print(eigen(V)$values)
      #   stop("V contains negative eigenvalues.")
      # }
      
      cholinv <- solve(t(chol(V)))
      Y2 <- matrix(cholinv%*%Y, ncol=1)
      X2 <- cholinv%*%X
      fit <- lm.fit(X2, Y2)
      beta.i <- matrix(fit$coefficients, ncol=1)
    }
    
    
    ## Defensive & debug conditions
    if(any(is.na(beta.i[c(which.fixed.cov, which.random.cov),]))) {
      print(beta.i)
      warning("For hl = ", hl," and vy = ", vy," the gls estimate of beta contains \"NA\". Consider using different hl or vy.")
      print(as.numeric(round(cbind(hl, vy, NA, t(beta1)), 4)))
      if (gridsearch){
        return(list(support = NA,
                    hl_vy = hl_vy))
      }else{
        return(list(support = NA,
                    hl_vy = hl_vy))
      }
    }
    
    ## Check for convergence
    con.count <- con.count + 1
    if (test.conv(beta.i, beta1, con.count, convergence)) {
      beta1<-beta.i
      break
    }else{
      beta1<-beta.i
    }
  }
  
  ## Compute residuals
  eY <- X%*%beta1
  resid1 <- Y-eY
  V.inverse <- solve(V)
  
  log.det.V <- mk.log.det.V(V = V, N = N)
  
  ## Log-likelihood
  sup1 <- -N/2*log(2*pi)-0.5*log.det.V-0.5*(t(resid1) %*% V.inverse%*%resid1)
  print(as.numeric(round(cbind(hl, vy, sup1, t(beta1)), 4)))
  
  if(gridsearch){
    return(list(support = sup1,
                hl_vy = hl_vy
                #V = V # This will cause memory overflow when N or Grid is large, or when less memory is available. O(N^2) + O(grid.vy * grid.hl)
                ))
  }else{
    beta1.var <- solve(t(X)%*%V.inverse%*%X)
    
    
    if(!is.null(fixed.cov) | !is.null(random.cov)){
      # ## Bias correction
      Vu <- diag(c(rep(0, N*length(beta1[-c(which.fixed.cov, which.random.cov),])), c(na.exclude(me.fixed.pred), na.exclude(me.pred))))
      Vd <- diag(c(rep(0, N*length(beta1[-c(which.fixed.cov, which.random.cov),])), c(sapply(Vd, function(e) diag(e)))))
      
      ## Center each column in X on its respective mean
      X0_intercept <- apply(matrix(X[, -c(which.fixed.cov, which.random.cov)], nrow=N), 2, function(e) e - mean(e))
      
      if(!is.null(fixed.cov)){
        X0_fixed <- apply(matrix(X[,which.fixed.cov], nrow=N), 2, function(e) e - mean(e))
      }else{
        X0_fixed <- NULL
      }
      
      if(!is.null(random.cov)){
        X0_random <- matrix(s.X, ncol=length(s.X), nrow=N)
      }else{
        X0_random <- NULL
      }
      
      X0 <- cbind(X0_intercept, X0_fixed, X0_random)
      
      correction <- matrix(Vu%*%pseudoinverse(Vd+Vu)%*%c(X0),  ncol=ncol(X), nrow=nrow(X), byrow=F)
      bias_corr <- pseudoinverse(t(X)%*%V.inverse%*%X)%*%t(X)%*%V.inverse%*%correction
      m<-length(beta1)
      K <- solve(diag(1,m,m)-bias_corr)
      
      beta1_bias_corr <- K%*%beta1
      beta1_bias_corr_var <- solve(K)%*%beta1.var%*%t(pseudoinverse(K))
    }

    
    ## 
    if(!is.null(random.cov)){
      X.ev <- calc.X(a = a, hl = hl, treepar, modelpar, seed, is.opt.reg = FALSE)
      ev.beta1.var <- pseudoinverse(t(X.ev)%*%V.inverse%*%X.ev)
      ev.beta1 <- ev.beta1.var%*%(t(X.ev)%*%V.inverse%*%Y)
      ev.reg <- list(coefficients = matrix(cbind(ev.beta1, sqrt(diag(ev.beta1.var))), nrow=ncol(X.ev), dimnames = list(colnames(X.ev), c("Estimates", "Std. error"))),
                     X = X.ev,
                     residuals = Y - (X.ev %*% ev.beta1))
    }else{
      ev.beta1.var <- NULL
      ev.beta1 <- NULL
      ev.reg <- NULL
    }
    opt.reg <- list(coefficients = matrix(cbind(beta1, sqrt(diag(beta1.var))), nrow=ncol(X), dimnames = list(colnames(X), c("Estimates", "Std. error"))),
                    X = X,
                    residuals = Y - (X %*% beta1))
    
    if(!is.null(fixed.cov) | !is.null(random.cov)){
      opt.reg$coefficients_bias_corr = matrix(cbind(beta1_bias_corr, sqrt(diag(beta1_bias_corr_var))), nrow=ncol(X), dimnames = list(colnames(X), c("Bias-corrected estimates", "Std. error")))
      opt.reg$residuals_bias_corr = Y - (X %*% beta1_bias_corr)
    }else{
      opt.reg$coefficients_bias_corr <- NULL
      opt.reg$residuals_bias_corr <- NULL
    }
    
    pred.mean <- X%*%beta1
    g.mean <- (t(rep(1, times=N))%*%solve(V)%*%Y)/sum(solve(V))
    sst<-t(Y-g.mean)%*% solve(V)%*%(Y-g.mean)
    sse<-t(Y-pred.mean)%*%solve(V)%*%(Y-pred.mean)
    r.squared<-(sst-sse)/sst
    
    return(list(support = sup1,
                V = V, 
                opt.reg = opt.reg,
                ev.reg = ev.reg,
                hl_vy = hl_vy,
                sst = sst,
                sse = sse,
                r.squared = r.squared))
  }
}