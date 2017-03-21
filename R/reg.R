## Generalized inverse, or pseudoinverse
pseudoinverse <-
  function (m){
    msvd <- svd(m)
    if (length(msvd$d) == 0) {
      return(array(0, dim(m)[2:1]))
    }
    else {
      return(msvd$v %*% (1 / msvd$d * t(msvd$u)))
    }
  }

## Design matrix
slouch.modelmatrix <- function(a, hl, tree, pars, control, seed, is.opt.reg = TRUE){
  list2env(tree, envir = environment())
  list2env(pars, envir = environment())
  list2env(control, envir = environment())
  
  list2env(seed, envir = environment())
  
  n <- length(tree$phy$tip.label)
  
  if(is.opt.reg == TRUE){
    rho <- (1-(1-exp(-a*T.term))/(a*T.term))
  }else{
    rho <- 1
  }
  
  if(!is.null(fixed.cov) | !is.null(random.cov)){
    covariates <- matrix(cbind(fixed.pred, rho*pred), nrow = n, dimnames = list(NULL, c(names.fixed.cov, names.random.cov)))
  }else{
    covariates <- NULL
  }

  if(estimate.bXa){
    if (!is.null(pars$random.cov) | !is.null(pars$fixed.cov)){
      bXa <- c(bXa = 1-exp(-a*T.term) - (1-(1-exp(-a*T.term))/(a*T.term)))
    }else{
      stop("bXa can not be estimated without continuous covariates.")
    }
  }else{
    bXa <- NULL
  }
  
  if(!is.null(fixed.fact)){
    w_regimes <- weight.matrix(tree$phy, a, lineages)
    X <- cbind(w_regimes, bXa, covariates)
  }else{
    if(!control$estimate.Ya){
      K <- cbind(Intercept = rep(1, length(tree$phy$tip.label)), 
                 bXa)
    }else{
      K <- cbind(Ya = exp(-a*T.term),
                 b0 = 1-exp(-a*T.term),
                 bXa)
    }
    X <- cbind(K, covariates)
  }
  return(X)
}

# Part "two" of variance-covariance matrix, hansen et al. 2008.
# It looks a little cryptic since it was vectorized from a nested for-loop
# Still is a bottleneck in performance. 01 sept 2016. Consider rewrite in Rcpp ?
calc.cm2 <- function(a, T.term, n, tia, tja, ta){
  ti <- matrix(rep(T.term, n), nrow=n)
  #tj <- t(ti)
  term0 <- exp(-a * tja) * (1 - exp(-a * ti)) / (a * ti)
  term2 <- (1 - exp(-a * ti))/(a * ti)
  num.prob <- ifelse(ta == 0, 1, (1 - exp(-a * ta)) / (a * ta))

  # res_old <- ((1-exp(-a*ti))/(a*ti))*(t((1-exp(-a*ti))/(a*ti))) -
  #   (exp(-a*tja)*(1-exp(-a*ti))/(a*ti) + t(exp(-a*tja)*(1-exp(-a*ti))/(a*ti)))*num.prob
  return(term2 * t(term2) - num.prob * (term0 + t(term0)))
}

varcov_model <- function(hl, vy, a, cm2, beta1, which.fixed.cov, which.random.cov, random.cov, T.term, fixed.cov, Vu_given_x, me.cov, n.pred, mecov.fixed.cov, n.fixed.pred, n, sigma_squared, ta, tij, me.response){
  
  ## Piece together V
  if (hl == 0){
    Vt <- diag(rep(vy, times = n))
  }else{
    if (!is.null(random.cov)){
      s1 <- sum(sigma_squared %*% ((beta1[which.random.cov, ])^2))
      Vt <- (s1 / (2 * a) + vy) * (1 - exp( -2 * a * ta)) * exp(-a * tij) + (s1 * ta * cm2)
    }else{
      Vt <- vy * (1 - exp(-2 * a * ta)) * exp(-a * tij)
    }
  }
  #print(paste0("alpha = ", a, "vy = ", vy, "beta1 = ", c(beta1)))
  return(Vt)
}

varcov_measurement <- function(tree, pars, control, seed, beta1, hl, a, which.fixed.cov, which.random.cov){
  list2env(tree, envir = environment())
  list2env(pars, envir = environment())
  list2env(control, envir = environment())
  
  list2env(seed, envir = environment())
  
  if (hl == 0 | is.null(random.cov)){
    rho2 <- 1
  }   else {
    rho2 <- (1 - (1 - exp(-a * T.term))/(a * T.term))^2
  }
  
  if(!is.null(fixed.cov) | !is.null(random.cov)){
    ## Update measurement error in X.
    ## Measurement error in predictor - take 2]
    beta2_Vu_given_x_list <-  mapply(function(Vu_given_xi, b, rho_squared) {Vu_given_xi * b^2 * rho_squared} , 
                                     Vu_given_xi = Vu_given_x, 
                                     b = beta1[c(which.fixed.cov, which.random.cov),], 
                                     rho_squared = c(lapply(seq_along(beta1[which.fixed.cov, ]), function(e) 1),
                                       lapply(seq_along(beta1[which.random.cov, ]), function(e) rho2)),
                                     SIMPLIFY = FALSE)
    beta2_Vu_given_x <- Reduce("+", beta2_Vu_given_x_list)
    
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
    beta2_Vu_given_x <- 0
    mcov <- 0
    mcov.fixed <- 0
  }
  
  V_me <- na.exclude(me.response) + beta2_Vu_given_x - mcov - mcov.fixed
  return(V_me)
}

reg <- function(hl_vy, tree, pars, control, seed, gridsearch = TRUE){
  list2env(tree, envir = environment())
  list2env(pars, envir = environment())
  list2env(control, envir = environment())
  list2env(seed, envir = environment())
  
  hl <- hl_vy[1]; vy <- hl_vy[2]
  n <- length(tree$phy$tip.label)
  
  if(hl == 0){
    a <- Inf ## When response is modeled only by fixed factors, of which one or more levels contain only internal regimes, Var(beta) becomes singular, throws error.
  }else{
    a <- log(2)/hl
    if (!is.null(random.cov)){
      cm2 <- calc.cm2(a, T.term, n, tia, tja, ta)
    }else cm2 <- NULL
  }
  X <- slouch.modelmatrix(a, hl, tree, pars, control, seed, is.opt.reg = TRUE)

  ## Initial OLS estimate of beta, disregarding variance covariance matrix
  beta1 <- matrix(lm.fit(X, Y)$coefficients, ncol=1)
  beta1.descriptor <- c(rep("Intercept", length(beta1) - length(names.fixed.cov) - length(names.random.cov)),
                        rep("Instantaneous cov", length(names.fixed.cov)),
                        rep("Random cov", length(names.random.cov)))
  which.fixed.cov <- which(beta1.descriptor == "Instantaneous cov")
  which.random.cov <- which(beta1.descriptor == "Random cov")

  con.count <- 0
  repeat{
    Vt <- varcov_model(hl, vy, a, cm2, beta1, which.fixed.cov, which.random.cov, random.cov, T.term, fixed.cov, Vu_given_x, me.cov, n.pred, mecov.fixed.cov, n.fixed.pred, n, sigma_squared, ta, tij, me.response)
    V_me <- varcov_measurement(tree, pars, control, seed, beta1, hl, a, which.fixed.cov, which.random.cov)
    V <- Vt + V_me

      # Linear transformation of both Y and model matrix, by the inverse of cholesky decomposition of V
      # Is this method ok to use? Seems to avoid some bugs (see below). R's lm.fit() does not throw
      # error when X is singular, but returns NA in the coefficients.
      # Problem with this transform is, does not work when matrix V is not positive semi-definite (Should this be possible?)
      #
      #  Error in solve.default(t(X) %*% V.inverse %*% X) : 
      #  system is computationally singular: reciprocal condition number = 2.28389e-22
      
      L <- t(chol(V))
      fit <- lm.fit(backsolve(L, X, upper.tri = FALSE), 
                    matrix(backsolve(L, Y, upper.tri = FALSE), ncol=1))
      beta.i <- matrix(fit$coefficients, ncol=1)
      
    ## Defensive & debug conditions
    if(any(is.na(beta.i[c(which.fixed.cov, which.random.cov),]))) {
      print(beta.i)
      warning("For hl = ", hl," and vy = ", vy," the gls estimate of beta contains \"NA\". Consider using different hl or vy.")
      print(as.numeric(round(cbind(hl, vy, NA, t(beta1)), 4)))
        return(list(support = NA,
                    hl_vy = hl_vy))
    }
      
    ## Check for convergence
    con.count <- con.count + 1
    if(all(abs(beta1[!is.na(beta1)] - beta.i[!is.na(beta.i)]) < convergence)) {
      beta1<-beta.i
      break
    }else{
      beta1<-beta.i
      
      if(con.count >= 50)
      {
        message("Warning, estimates did not converge after 50 iterations, last estimates printed out")
        break
      }
    }
  }
  
  log.det.V = 2*sum(log(diag(L))) ## Special case for positive definite matrix V, more numerically stable than (log*det(V)) routine for large V
  
  ## Log-likelihood
  sup1 <- - n/2*log(2*pi) - 0.5*log.det.V - 0.5*crossprod(fit$residuals)
  if(verbose){
    print(as.numeric(round(cbind(hl, vy, sup1, t(beta1)), 4)))
  }
  
  if(gridsearch){
    return(list(support = sup1,
                hl_vy = hl_vy
                #V = V # This will cause memory overflow when n or Grid is large, or when less memory is available. O(n^2) + O(grid.vy * grid.hl)
                ))
  }else{
    
    ## Following is less optimized, more identical to publication syntax because it only has to run once.
    V.inverse <- solve(V)
    beta1.var <- solve(t(X)%*%V.inverse%*%X)
    
    
    if(!is.null(fixed.cov) | !is.null(random.cov)){
      # ## Bias correction
      Vu <- diag(c(rep(0, n*length(beta1[-c(which.fixed.cov, which.random.cov),])), c(na.exclude(me.fixed.pred), na.exclude(me.pred))))
      Vd <- diag(c(rep(0, n*length(beta1[-c(which.fixed.cov, which.random.cov),])), c(sapply(Vd, function(e) diag(e)))))
      
      ## Center each column in X on its respective mean
      X0_intercept <- apply(matrix(X[, -c(which.fixed.cov, which.random.cov)], nrow=n), 2, function(e) e - mean(e))
      
      if(!is.null(fixed.cov)){
        X0_fixed <- apply(matrix(X[,which.fixed.cov], nrow=n), 2, function(e) e - mean(e))
      }else{
        X0_fixed <- NULL
      }
      
      if(!is.null(random.cov)){
        X0_random <- matrix(sigma_squared, ncol=length(sigma_squared), nrow=n, byrow = TRUE)
      }else{
        X0_random <- NULL
      }
      
      X0 <- cbind(X0_intercept, X0_fixed, X0_random)
      
      correction <- matrix(Vu%*%pseudoinverse(Vd+Vu)%*%c(X0),  ncol=ncol(X), nrow=nrow(X), byrow=F)
      bias_corr <- pseudoinverse(t(X)%*%V.inverse%*%X)%*%t(X)%*%V.inverse%*%correction
      m<-length(beta1)
      K <- solve(diag(1,m,m)-bias_corr)
      
      beta1_bias_corr <- K%*%beta1
      #beta1_bias_corr_var <- solve(K)%*%beta1.var%*%t(pseudoinverse(K))
      beta1_bias_corr_var <- solve(K)%*%beta1.var%*%t(solve(K))
    }

    
    ## 
    if(!is.null(random.cov)){
      X.ev <- slouch.modelmatrix(a = a, hl = hl, tree, pars, control, seed, is.opt.reg = FALSE)
      ev.beta1.var <- pseudoinverse(t(X.ev)%*%V.inverse%*%X.ev)
      ev.beta1 <- ev.beta1.var%*%(t(X.ev)%*%V.inverse%*%Y)
      ev.reg <- list(coefficients = matrix(cbind(ev.beta1, sqrt(diag(ev.beta1.var))), 
                                           nrow=ncol(X.ev), 
                                           dimnames = list(colnames(X.ev), c("Estimates", "Std. error"))),
                     X = X.ev,
                     residuals = Y - (X.ev %*% ev.beta1))
    }else{
      ev.beta1.var <- NULL
      ev.beta1 <- NULL
      ev.reg <- NULL
    }
    opt.reg <- list(coefficients = matrix(cbind(beta1, sqrt(diag(beta1.var))), 
                                          nrow=ncol(X), 
                                          dimnames = list(colnames(X), c("Estimates", "Std. error"))),
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
    g.mean <- (t(rep(1, times=n))%*%solve(V)%*%Y)/sum(solve(V))
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