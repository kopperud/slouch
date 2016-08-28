calc.X <- function(a, hl, treepar, modelpar, seed, is.opt.reg){
  list2env(modelpar, envir = environment())
  list2env(treepar, envir = environment())
  list2env(seed, envir = environment())
  
  if(is.opt.reg == TRUE){
    rho <- (1-(1-exp(-a*T.term))/(a*T.term))
  }else{
    rho <- 1
  }
  if(!is.null(fixed.fact)){
    cbind(weight.matrix(a, topology, times, N, regime.specs, fixed.pred, intercept), rho*pred)
  }else{
    if(!is.null(modelpar$intercept) | (is.null(modelpar$random.cov) & is.null(modelpar$fixed.cov))){
      K <- 1
    }else{
      K <- cbind(exp(-a*T.term),
                 1-exp(-a*T.term),
                 if (!is.null(modelpar$random.cov) | !is.null(modelpar$fixed.cov)) 1-exp(-a*T.term)-(1-(1-exp(-a*T.term))/(a*T.term)) else NULL
      )
    }
    matrix(cbind(K,
                 fixed.pred,
                 rho*pred), 
           nrow=N)
  }
}


regression.closures <- function(treepar, modelpar, seed){
  ## bind global variables to this environment
  list2env(treepar, envir = environment())
  list2env(modelpar, envir = environment())
  list2env(seed, envir = environment())
  
  ## Function to calculate covariances between response and the stochastic predictor, to be subtracted in the diagonal of V
  ## BROKEN! Same as above. Needs fix.
  calc.mcov <- function(a, beta1){
    if(sum(me.cov) == 0){
      0
    }else{
      diag(rowSums(matrix(data=as.numeric(me.cov)*t(kronecker(2*beta1[which.random.cov,],(1-(1-exp(-a*T.term))/(a*T.term)))), ncol=n.pred)))
    }
  }
  
  ## Function to calculate covariances between response and instantaneous predictor, to be subtracted in the diagonal of V
  calc.mcov.fixed <- function(a, beta1){
    if(sum(mecov.fixed.cov) == 0){
      matrix(0, nrow=N, ncol=N)
    }else{
      diag(rowSums(matrix(data=as.numeric(mecov.fixed.cov)*t(kronecker(2*beta1[which.fixed.cov,],(1-(1-exp(-a*T.term))/(a*T.term)))), ncol=n.fixed.pred)))
    }
  }
  
  ## Function to calculate V
  calc.V <- function(hl, vy, a, cm2, beta1){
    
    ## Update measurement error in X
    if (hl == 0 | is.null(random.cov)){
      y <- 1
    }   else {
      y <- ((1-(1-exp(-a*T.term))/(a*T.term))*(1-(1-exp(-a*T.term))/(a*T.term)))
    }
    
    ## measurement error of predictor variables
    if(!is.null(fixed.cov) | !is.null(random.cov)){
      obs_var_con <- matrix(0, nrow=N, ncol=N)
      for (e in seq(from=1, to=ncol(x.ols), by=1)){
        for (j in seq(from=1, to=ncol(x.ols), by=1)) {
          tmp <- error_condition[xx[e]:(e*N),xx[j]:(j*N)]*beta1[e]*beta1[j]*y
          obs_var_con <- obs_var_con + tmp
          

        }
      }
    mepredictorf <- function(){
      obs_var_con <- matrix(0, nrow=N, ncol=N)
      for (e in seq(from=1, to=ncol(x.ols), by=1)){
        for (j in seq(from=1, to=ncol(x.ols), by=1)) {
          tmp <- error_condition[xx[e]:(e*N),xx[j]:(j*N)]*beta1[e]*beta1[j]*y
          obs_var_con <- obs_var_con + tmp
        }
      }
    }

    print(microbenchmark(mepredictorf()))
      
      
      # Covariances
      mcov <- calc.mcov(a, beta1)
      mcov.fixed <- calc.mcov.fixed(a, beta1)
    }else{
      obs_var_con <- 0
      mcov <- 0
      mcov.fixed <- 0
    }
    
    ## Piece together V
    if(hl == 0){
      cm0 <- diag(rep(vy, times=N))
    }else{
      if(!is.null(random.cov)){
        s1 <- as.numeric(s.X%*%((beta1[which.random.cov,])^2))
        cm0 <- (s1/(2*a)+vy)*(1-exp(-2*a*ta))*exp(-a*tij) + (s1*ta*cm2)
      }else{
        cm0 <- vy*(1-exp(-2*a*ta))*exp(-a*tij)
      }
    }
    V <- cm0 + na.exclude(me.response) + obs_var_con - mcov - mcov.fixed
    return(V)
  }
  
  # Part "two" of variance-covariance matrix, hansen et al. 2008. Everything after sigma^2_theta * ta ????? Looks like equation is modified to 
  # account for non-ultrametric trees.
  # It looks a little cryptic since it was vectorized from a nested for-loop
  calc.cm2 <- function(a){
    T.row <- replicate(N,T.term)
    T.col <- t(T.row)
    num.prob <- ifelse(ta == 0, 1, (1-exp(-a*ta))/(a*ta))
    return(((1-exp(-a*T.row))/(a*T.row))*((1-exp(-a*T.col))/(a*T.col)) - (exp(-a*tia)*(1-exp(-a*T.row))/(a*T.col) + exp(-a*tja)*(1-exp(-a*T.row))/(a*T.row))*num.prob)
  }
  
  ## General test for beta convergence
  test.conv <- function(beta.i, beta1, con.count){
    if(!is.null(modelpar$intercept) | (is.null(fixed.cov) & is.null(random.cov))){
      test <- ifelse(abs(as.numeric(beta.i - beta1)) <= convergence, 0, 1)
    }else{
      #test <- ifelse(abs(as.numeric(beta.i - beta1))[-(1:2)] <= convergence, 0, 1)
      test <- ifelse(abs(as.numeric(beta.i - beta1)) <= convergence, 0, 1)
      ## Effectively removes beta[1:2] <= 0.001 from being criteria in convergence, when non-ultrametric with continuous covariates
    }
    
    if(sum(test)==0) return (TRUE)
    if(con.count >= 50)
    {
      message("Warning, estimates did not converge after 50 iterations, last estimates printed out")
      return(TRUE)
    }
    return(FALSE)
  }
  
  ## Function for regression & grid search
  slouch.regression <- function(hl_vy){
    hl <- hl_vy[1]; vy <- hl_vy[2]
    
    if(hl == 0){
      a <- Inf ## When response is modeled only by fixed factors, of which one or more levels contain only internal regimes, Var(beta) becomes singular, throws error.
    }else{
      a <- log(2)/hl
      if (!is.null(random.cov)){
        cm2 <- calc.cm2(a)
        #print(microbenchmark(calc.cm2(a)))
      }else cm2 <- NULL
    }
    #X <- calc.X(a, hl)
    X <- calc.X(a, hl, treepar, modelpar, seed, is.opt.reg = TRUE)
    #print(microbenchmark(calc.X(a, hl, treepar, modelpar, seed, is.opt.reg = TRUE)))
    
    ## ols.beta1 exists from OLS estimate
    beta1 <- ols.beta1
    
    con.count <- 0
    repeat{
      V <- calc.V(hl, vy, a, cm2, beta1 = beta1)
      #print(microbenchmark(calc.V(hl, vy, a, cm2, beta1 = beta1)))
      
      V.inverse<-solve(V)
      #beta.i.var <- solve(t(X)%*%V.inverse%*%X)
      beta.i.var <- pseudoinverse(t(X)%*%V.inverse%*%X)
      
      #### Ask Thomas about this one.
      if(Inf %in% beta.i.var) {print("Pseudoinverse of (XT * V * X) contained values = Inf, which were set to 10^300")}
      beta.i.var <-replace(beta.i.var, beta.i.var ==Inf, 10^300)
      if(-Inf %in% beta.i.var) {print("Pseudoinverse of (XT * V * X) contained values = -Inf, which were set to -10^300")}
      beta.i.var <-replace(beta.i.var, beta.i.var ==-Inf, -10^300)
      
      beta.i<-beta.i.var%*%(t(X)%*%V.inverse%*%Y)
      #print(beta.i);print(beta1)
      
      ## Check for convergence
      con.count <- con.count + 1
      if (test.conv(beta.i, beta1, con.count)) {
        beta1<-beta.i
        break
      }else{
        beta1<-beta.i
      }
    }
    
    ## Compute residuals
    eY<-X%*%beta1
    resid1<-Y-eY
    
    ## Log-likelihood
    log.det.V <- mk.log.det.V(V = V, N = N)
    #print(microbenchmark(mk.log.det.V(V = V, N = N)))
    sup1 <- -N/2*log(2*pi)-0.5*log.det.V-0.5*(t(resid1) %*% V.inverse%*%resid1)
    print(as.numeric(round(cbind(hl, vy, sup1, t(beta1)), 4)))
    # printmsg <- as.numeric(round(cbind(hl, vy, sup1, t(beta1)), 4))
    # cat(printmsg); cat("\n")
    
    
    list(support = sup1,
         V = V, # This will cause memory overflow when N or Grid is large, or when less memory is available. O(N^2) + O(grid.vy * grid.hl)
         beta1 = beta1,
         X = X,
         Y = Y,
         beta1.var = beta.i.var,
         alpha.est = a,
         vy.est = vy)
  }
  
  all.closures <- list(#calc.X = calc.X,
                       calc.mcov = calc.mcov,
                       calc.mcov.fixed = calc.mcov.fixed,
                       calc.V = calc.V,
                       slouch.regression = slouch.regression)
  return(all.closures)
}

