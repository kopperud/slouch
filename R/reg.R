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
slouch.modelmatrix <- function(a, hl, tree, observations, control, is.opt.reg = TRUE){

  
  n <- length(tree$phy$tip.label)
  
  if(is.opt.reg == TRUE){
    rho <- (1 - (1 - exp(-a * tree$T.term))/(a * tree$T.term))
  }else{
    rho <- 1
  }
  
  if(!is.null(observations$fixed.cov) | !is.null(observations$random.cov)){
    covariates <- matrix(cbind(observations$fixed.cov, rho * observations$random.cov), 
                         nrow = n, 
                         dimnames = list(NULL, c(observations$names.fixed.cov, 
                                                 observations$names.random.cov)))
  }else{
    covariates <- NULL
  }

  if(control$estimate.bXa){
    if (!is.null(observations$random.cov) | !is.null(observations$fixed.cov)){
      bXa <- c(bXa = 1 - exp(-a * tree$T.term) - (1 - (1 - exp(-a * tree$T.term))/(a * tree$T.term)))
    }else{
      stop("bXa can not be estimated without continuous covariates.")
    }
  }else{
    bXa <- NULL
  }
  
  if(!is.null(observations$fixed.fact)){
    w_regimes <- weight.matrix(tree$phy, a, tree$lineages)
    X <- cbind(w_regimes, bXa, covariates)
  }else{
    if(!control$estimate.Ya){
      K <- cbind(Intercept = rep(1, length(tree$phy$tip.label)), 
                 bXa)
    }else{
      K <- cbind(Ya = exp(-a * tree$T.term),
                 b0 = 1-exp(-a * tree$T.term),
                 bXa)
    }
    X <- cbind(K, covariates)
  }
  return(X)
}

#varcov_model <- function(hl, vy, a, beta1, which.fixed.cov, which.random.cov, random.cov, T.term, fixed.cov, Vu_given_x, mecov.random.cov, mecov.fixed.cov, n, sigma_squared, phy, ta, tij, tja, me.response){
varcov_model <- function(hl, sigma2_y, a, beta1, which.fixed.cov, which.random.cov, tree, observations, seed){

  tij <- tree$tij
  tja <- tree$tja
  ta <- tree$ta
  n <- length(tree$phy$tip.label)
  
  if (hl == 0){
    Vt <- diag(rep(sigma2_y, times = n))
    #print(Vt)
  }else if (a < 1e-14){
    Vt <- sigma2_y * ta
  }else{
    if (!is.null(observations$random.cov)){
      s1 <- sum(seed$sigma_squared %*% ((beta1[which.random.cov, ])^2))
      ti <- matrix(rep(tree$T.term, n), nrow=length(tree$phy$tip.label))
      
      term0 <- ((s1 + sigma2_y) / (2 * a)) * (1 - exp( -2 * a * ta)) * exp(-a * tij)
      term1 <- (1 - exp(-a * ti))/(a * ti) 
      term2 <- exp(-a * tja) * (1 - exp(-a * ti)) / (a * ti) 
      
      Vt <- term0 + s1*(ta * term1 * t(term1) - ((1 - exp(-a * ta))/a) * (term2 + t(term2)))
      
    }else{
      #Vt <- vy * (1 - exp(-2 * a * ta)) * exp(-a * tij)
      Vt <- sigma2_y/(2*a) * observations$closures$V_fixed_partial(a)
    }
  }
  return(Vt)
}

varcov_measurement <- function(observations, seed, beta1, hl, a, T.term, which.fixed.cov, which.random.cov){


  rho2 <- (1 - (1 - exp(-a * T.term))/(a * T.term))^2
  
  if(length(which.fixed.cov) > 0 | length(which.random.cov) > 0){
    ## Update measurement error in X.
    ## Measurement error in predictor - take 2]
    beta2_Vu_given_x_list <-  mapply(function(Vu_given_xi, b, rho_squared) {Vu_given_xi * b^2 * rho_squared} , 
                                     Vu_given_xi = seed$Vu_given_x, 
                                     b = beta1[c(which.fixed.cov, which.random.cov),], 
                                     rho_squared = c(lapply(seq_along(beta1[which.fixed.cov, ]), function(e) 1),
                                                     lapply(seq_along(beta1[which.random.cov, ]), function(e) rho2)),
                                     SIMPLIFY = FALSE)
    beta2_Vu_given_x <- Reduce("+", beta2_Vu_given_x_list)
    
    # ## BROKEN!? Needs test.
    # calculate covariances between response and the stochastic predictor, to be subtracted in the diagonal of V
    if(sum(observations$mecov.random.cov) == 0){
      mcov <- 0
    }else{
      mcov <- rowSums(matrix(data=as.numeric(observations$mecov.random.cov)*t(kronecker(2*beta1[which.random.cov,],(1-(1-exp(-a*T.term))/(a*T.term)))), ncol = ncol(observations$random.cov)))
    }
    
    if(sum(observations$mecov.fixed.cov) == 0){
      mcov.fixed <- 0
    }else{
      mcov.fixed <- rowSums(matrix(data=as.numeric(observations$mecov.fixed.cov)*t(kronecker(2*beta1[which.fixed.cov,],(1-(1-exp(-a*T.term))/(a*T.term)))), ncol =  ncol(observations$fixed.cov)))
    }
  }else{
    beta2_Vu_given_x <- 0
    mcov <- 0
    mcov.fixed <- 0
  }
  
  V_me <- observations$me.response + beta2_Vu_given_x - mcov - mcov.fixed
  return(V_me)
}

reg <- function(par, tree, observations, control, seed, gridsearch = TRUE){
  Y <- observations$Y

  if(is.numeric(par)){
    if(length(par) >= 2){
      par <- as.list(par)
    }else{
      par <- list(sigma2_y = par)
    }
  }
  
  if (control$model == "ou"){
    hl <- par$hl

    if(hl == 0){
      a <- Inf ## When response is modeled only by fixed factors, of which one or more levels contain only internal regimes, Var(beta) becomes singular, throws error.
    }else{
      a <- log(2)/hl
    }
    
    if(!is.null(par$vy)){
      vy <- par$vy
      sigma2_y <- vy*2*a 
    }else{
      sigma2_y <- par$sigma2_y
      vy <- sigma2_y/(2*a)
    }
    

  }else{
    #print(par)
    #sigma2_y <- if(length(par) == 1) par[[1]] else par$sigma2_y
    sigma2_y <- par$sigma2_y
    #sigma2_y <- unlist(par)
    hl <- Inf
    a <- 0
    vy <- NULL
  }
  

  
  
  n <- length(tree$phy$tip.label)
  

  
  if (control$verbose){
    cat(as.numeric(round(cbind(hl, (if(!is.null(vy)) vy else sigma2_y)), 4)))
    cat(" ")
  }
  
  X <- slouch.modelmatrix(a, hl, tree, observations, control, is.opt.reg = TRUE)


  ## Initial OLS estimate of beta, disregarding variance covariance matrix
  beta1 <- matrix(stats::lm.fit(X, Y)$coefficients, ncol=1)
  beta1.descriptor <- c(rep("Intercept", length(beta1) - length(observations$names.fixed.cov) - length(observations$names.random.cov)),
                        rep("Instantaneous cov", length(observations$names.fixed.cov)),
                        rep("Random cov", length(observations$names.random.cov)))
  which.fixed.cov <- which(beta1.descriptor == "Instantaneous cov")
  which.random.cov <- which(beta1.descriptor == "Random cov")

  
  con.count <- 0
  repeat{
    Vt <- varcov_model(hl, sigma2_y, a, beta1, which.fixed.cov, which.random.cov, tree, observations, seed)
    V_me <- varcov_measurement(observations, seed, beta1, hl, a, tree$T.term, which.fixed.cov, which.random.cov)
    V <- Vt + diag(V_me)
    
    L <- t(chol(V))

    # Note: Calculation of regression coefficients deviates from paper notation, but should be algebraically 
    # equivalent. This is done for an increase in computational speed, convenience of diagnosing problems
    # with better error-messages, and it is possibly more numerically stable.
    #
    # Instead of writing out the naive GLS beta estimator as-is, we have
    #
    # Linear transformation of both Y and model matrix, by the inverse of cholesky decomposition of V
    # With this method, the program will not crash even if X is singular. In this case, R's lm.fit() does not throw
    # error, but return coefficients that contain "NA", and subsequently the log-likelihood is evaluated as "NA".

    fit <- lm.fit(backsolve(L, X, upper.tri = FALSE), 
                  matrix(backsolve(L, Y, upper.tri = FALSE), ncol=1))
    beta.i <- matrix(fit$coefficients, ncol=1)
    
    ## Defensive & debug conditions
    if(any(is.na(beta.i[c(which.fixed.cov, which.random.cov), ]))) {
      cat("\n")
      print(beta.i)
      warning("For hl = ", hl," and vy = ", vy," the gls estimate of beta contains \"NA\". Consider using different hl or vy.")
      print(as.numeric(round(cbind(hl, vy, NA, t(beta1)), 4)))
        return(list(support = NA,
                    par = par))
    }
      
    ## Check for convergence
    con.count <- con.count + 1
    if(all(abs(beta1[!is.na(beta1)] - beta.i[!is.na(beta.i)]) < control$convergence)) {
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
  if(control$verbose){
    cat(as.numeric(round(cbind(sup1, t(beta1)), 4)))
    cat("\n")
  }
  
  
  if(gridsearch){
    return(list(support = sup1,
                par = par
                #V = V # This will cause memory overflow when n or Grid is large, or when less memory is available. O(n^2) + O(grid.vy * grid.hl)
                ))
  }else{
    
    ## Following is less optimized, more identical to publication syntax because it only has to run once.
    V.inverse <- solve(V)
    beta1.var <- solve(t(X)%*%V.inverse%*%X)
    

    
    ## Evolutionary regression
    X0 <- slouch.modelmatrix(a = a, hl = hl, tree, observations, control, is.opt.reg = FALSE)
    if(!is.null(observations$random.cov)){
      ev.beta1.var <- pseudoinverse(t(X0)%*%V.inverse%*%X0)
      ev.beta1 <- ev.beta1.var%*%(t(X0)%*%V.inverse%*%Y)
      ev.reg <- c(list(coefficients = matrix(cbind(ev.beta1, sqrt(diag(ev.beta1.var))), 
                                             nrow=ncol(X0), 
                                             dimnames = list(colnames(X0), c("Estimates", "Std. error"))),
                       X = X0,
                       residuals = Y - (X0 %*% ev.beta1)),
                  bias_correction(ev.beta1, ev.beta1.var, Y, X0, V, which.fixed.cov, which.random.cov, seed))
    }else{
      ev.beta1.var <- NULL
      ev.beta1 <- NULL
      ev.reg <- NULL
    }
    
    opt.reg <- c(list(coefficients = matrix(cbind(beta1, sqrt(diag(beta1.var))), 
                                            nrow=ncol(X), 
                                            dimnames = list(colnames(X), c("Estimates", "Std. error"))),
                      X = X,
                      residuals = Y - (X0 %*% beta1)),
                 if (!is.null(observations$fixed.cov) | !is.null(observations$random.cov)) bias_correction(beta1, beta1.var, Y, X, V, which.fixed.cov, which.random.cov, seed) else NULL)
    
    pred.mean <- X%*%beta1
    g.mean <- (t(rep(1, times = n)) %*% solve(V) %*% Y) / sum(solve(V))
    sst <- t(Y - c(g.mean)) %*% solve(V) %*% (Y - c(g.mean))
    sse <- t(Y - pred.mean) %*% solve(V) %*% (Y - pred.mean)
    r.squared <- (sst - sse) / sst
    
    return(list(support = sup1,
                par = par,
                V = V, 
                opt.reg = opt.reg,
                ev.reg = ev.reg,
                sst = sst,
                sse = sse,
                r.squared = r.squared))
  }
}

