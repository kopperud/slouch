regression.closures <- function(treepar, modelpar, seed){
  ## bind global variables to this environment
  list2env(treepar, envir = environment())
  list2env(modelpar, envir = environment())
  list2env(seed, envir = environment())
  
  
  
  ## Establish function to calculate design matrix X
  if(!is.null(fixed.fact)){
    calc.X <- function(a, hl){
      if(hl == 0){
        cbind(weight.matrix(a, topology, times, N, regime.specs, fixed.pred, intercept), pred)
      }else{
        cbind(weight.matrix(a, topology, times, N, regime.specs, fixed.pred, intercept), (1-(1-exp(-a*T.term))/(a*T.term))*pred)
      }
    }
  }else{
    calc.X <- function(a, hl){
      if(hl == 0){
        cbind(1, fixed.pred, pred)
      }else{
        if (ultrametric == TRUE){
          matrix(cbind(1, fixed.pred, (1-(1-exp(-a*T.term))/(a*T.term))*pred), nrow=N)
        }else{
          matrix(cbind(1-exp(-a*T.term), 
                       if (is.null(random.cov) & is.null(fixed.cov)) 1-exp(-a*T.term)-(1-(1-exp(-a*T.term))/(a*T.term)) else NULL,
                       exp(-a*T.term),
                       fixed.pred,
                       (1-(1-exp(-a*T.term))/(a*T.term))*pred),
                nrow=N)
        }
      }
    }
  }
  
  
  ## Function to calculate variance of instantaneous predictor in the OU-process
  calc.s1 <- function(beta1){
    if(ultrametric == TRUE){
      as.numeric(s.X%*%(beta1[(2+n.fixed.pred):(n.pred+n.fixed.pred+1),]*beta1[(2+n.fixed.pred):(n.pred+n.fixed.pred+1),]))
    }else{
      as.numeric(s.X%*%(beta1[(4+n.fixed.pred):(n.pred+n.fixed.pred+3),]*beta1[(4+n.fixed.pred):(n.pred+n.fixed.pred+3),]))
    }
  }
  
  ## Function to calculate covariances of the stochastic predictor, to be subtracted in the diagonal of V
  calc.mcov <- function(a, beta1){
    if(!is.null(random.cov)){
      if(ultrametric == TRUE){
        diag(rowSums(matrix(data=as.numeric(me.cov)*t(kronecker(2*beta1[(2+n.fixed.pred):(n.pred+n.fixed.pred+1),], (1-(1-exp(-a*T.term))/(a*T.term)))), ncol=n.pred)))
      }else{
        diag(rowSums(matrix(data=as.numeric(me.cov)*t(kronecker(2*beta1[(4+n.fixed.pred):(n.pred+n.fixed.pred+3),], (1-(1-exp(-a*T.term))/(a*T.term)))), ncol=n.pred)))
      }
    }else{
      0
    }
  }
  
  ## Function to calculate covariances of the instantaneous predictor, to be subtracted in the diagonal of V
  calc.mcov.fixed <- function(a, beta){
    if(sum(me.fixed.cov) == 0){
      matrix(0, nrow=N, ncol=N)
    }else{
      if(ultrametric == TRUE){
        diag(as.numeric(me.fixed.cov%*%(2*beta1[(2 + n.factor):(length(beta1)-n.pred),])))
      }else{
        diag(as.numeric(me.fixed.cov%*%(2*beta1[(4 + n.factor):(length(beta1)-n.pred),])))
      }
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
    
    if(!is.null(fixed.cov) | !is.null(random.cov)){
      obs_var_con <- matrix(0, nrow=N, ncol=N)
      for (e in seq(from=1, to=ncol(x.ols), by=1)){
        for (j in seq(from=1, to=ncol(x.ols), by=1)) {
          tmp <- error_condition[xx[e]:(e*N),xx[j]:(j*N)]*beta1[e]*beta1[j]*y
          obs_var_con <- obs_var_con + tmp
        }
      }
      
      mcov <- calc.mcov(a, beta1)
      mcov.fixed <- calc.mcov.fixed(a, beta1)
    }else{
      obs_var_con <- 0
      mcov <- 0
      mcov.fixed <- 0
    }
      
    ## Piece together V
    if(hl == 0){
      a <- Inf

      cm0 <- diag(rep(vy, times=N))
    }else{
      a <- log(2)/hl


      if(!is.null(random.cov)){
        s1 <- calc.s1(beta1)
        cm0 <- (s1/(2*a)+vy)*(1-exp(-2*a*ta))*exp(-a*tij) + (s1*ta*cm2)
      }else{
        cm0 <- vy*(1-exp(-2*a*ta))*exp(-a*tij)
      }
    }
    V <- cm0 + na.exclude(me.response) + obs_var_con - mcov - mcov.fixed
    return(V)
  }
  
  calc.cm2 <- function(a){
    T.row <- replicate(N,T.term)
    T.col <- t(T.row)
    num.prob <- ifelse(ta == 0, 1, (1-exp(-a*ta))/(a*ta))
    return(((1-exp(-a*T.row))/(a*T.row))*((1-exp(-a*T.col))/(a*T.col)) - (exp(-a*tia)*(1-exp(-a*T.row))/(a*T.col) + exp(-a*tja)*(1-exp(-a*T.row))/(a*T.row))*num.prob)
  }
  
  ## Function for regression & grid search
  slouch.regression <- function(hl_vy){
    hl <- hl_vy[1]; vy <- hl_vy[2]
    
    if(hl == 0){
      a <- Inf
    }else{
      a <- log(2)/hl
      if (!is.null(random.cov)){
        cm2 <- calc.cm2(a)
      }else cm2 <- NULL
    }
    X <- calc.X(a, hl)
    
    ## ols.beta1 exists from OLS estimate
    beta1 <- ols.beta1
    
    con.count <- 0
    repeat{
      V <- calc.V(hl, vy, a, cm2, beta1 = beta1) # beta1, cm2
      
      
      V.inverse<-solve(V)
      beta.i.var <- pseudoinverse(t(X)%*%V.inverse%*%X)
      beta.i<-beta.i.var%*%(t(X)%*%V.inverse%*%Y)
      
      ## Check for convergence
      con.count <- con.count + 1
      if (test.conv(beta.i, beta1, convergence, con.count, ultrametric)) break
      beta1<-beta.i
    }
    # Remember to bring the last beta
    beta1<-beta.i
    
    ## Compute residuals
    eY<-X%*%beta1
    resid1<-Y-eY
    
    log.det.V <- mk.log.det.V(V = V, N = N)
    
    sup1 <- -N/2*log(2*pi)-0.5*log.det.V-0.5*(t(resid1) %*% V.inverse%*%resid1)
    print(as.numeric(round(cbind(hl, vy, sup1, t(beta1)), 4)))
    
    list(support = sup1,
         V = V,
         beta1 = beta1,
         X = X,
         Y = Y,
         beta1.var = beta.i.var,
         alpha.est = a,
         vy.est = vy)
  }
  
  all.closures <- list(calc.X = calc.X,
                       calc.s1 = calc.s1,
                       calc.mcov = calc.mcov,
                       calc.mcov.fixed = calc.mcov.fixed,
                       calc.V = calc.V,
                       slouch.regression = slouch.regression)
  return(all.closures)
}
