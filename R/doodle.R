regression.closures <- function(treepar, modelpar, seed){
  ## bind global variables to this environment
  for (i in c(treepar, modelpar, seed)){
    list2env(i, envir = environment())
  }
  
  ## Establish function to calculate design matrix X
  if(is.null(fixed.fact) == FALSE){
    calc.X <- function(a){
      if(hl == 0){
        cbind(weight.matrix(a, topology, times, N, regime.specs, fixed.pred, intercept), pred)
      }else{
        cbind(weight.matrix(a, topology, times, N, regime.specs, fixed.pred, intercept), (1-(1-exp(-a*T.term))/(a*T.term))*pred)
      }
    }
  }else{
    calc.X <- function(a){
      if(hl == 0){
        X<-cbind(1, fixed.pred, pred)
      }else{
        if (ultrametric == TRUE){
          cbind(1, fixed.pred, (1-(1-exp(-a*T.term))/(a*T.term))*pred)
        }else{
          cbind(1-exp(-a*T.term), 1-exp(-a*T.term)-(1-(1-exp(-a*T.term))/(a*T.term)), exp(-a*T.term), fixed.pred, (1-(1-exp(-a*T.term))/(a*T.term))*pred)
        }
      }
    }
  }
  
  ## Function to calculate variance of predictor in the OU-process
  calc.s1 <- function(beta1){
    if(ultrametric == TRUE){
      as.numeric(s.X%*%(beta1[(2+n.fixed.pred):(n.pred+n.fixed.pred+1),]*beta1[(2+n.fixed.pred):(n.pred+n.fixed.pred+1),]))
    }else{
      as.numeric(s.X%*%(beta1[(4+n.fixed.pred):(n.pred+n.fixed.pred+3),]*beta1[(4+n.fixed.pred):(n.pred+n.fixed.pred+3),]))
    }
  }
  
  ## Function to calculate covariances of the stochastic predictor, to be subtracted in the diagonal of V
  calc.mcov <- function(a, beta1){
    if(ultrametric == TRUE){
      diag(rowSums(matrix(data=as.numeric(me.cov)*t(kronecker(2*beta1[(2+n.fixed.pred):(n.pred+n.fixed.pred+1),], (1-(1-exp(-a*T.term))/(a*T.term)))), ncol=n.pred)))
    }else{
      diag(rowSums(matrix(data=as.numeric(me.cov)*t(kronecker(2*beta1[(4+n.fixed.pred):(n.pred+n.fixed.pred+3),], (1-(1-exp(-a*T.term))/(a*T.term)))), ncol=n.pred)))
    }
  }
  
  ## Function to calculate covariances of the instantaneous predictor, to be subtracted in the diagonal of V
  calc.mcov.fixed <- function(a, beta){
    if(ultrametric == TRUE){
      diag(as.numeric(me.fixed.cov%*%(2*beta1[2:(length(beta1)-n.pred),])))
    }else{
      diag(as.numeric(me.fixed.cov%*%(2*beta1[4:(length(beta1)-n.pred),])))
    }
  }
  
  ## Function to calculate V
  calc.V <- function(hl, vy, beta1, cm2){
    
    ## Update measurement error in X
    obs_var_con <- matrix(0, nrow=N, ncol=N)
    for (e in seq(from=1, to=ncol(x.ols), by=1)){
      for (j in seq(from=1, to=ncol(x.ols), by=1)) {
        tmp <- error_condition[xx[e]:(e*N),xx[j]:(j*N)]*beta1[e]*beta1[j]
        obs_var_con <- obs_var_con + tmp
      }
    }
    
    ## Piece together V
    if(hl == 0){
      a <- Inf
      V<-diag(rep(vy, times=N))+na.exclude(me.response)+ obs_var_con - 
        diag(as.numeric(me.cov%*%(2*beta1[(n.fixed+1+n.fixed.pred):length(beta1),]))) -
        diag(as.numeric(me.fixed.cov%*%(2*beta1[(n.fixed+1):(length(beta1)-n.pred),])))
    }else{
      a <- log(2)/hl
      s1 <- make.s1(beta1)
      mcov <- make.mcov(a, beta1)
      mcov.fixed <- make.mcov.fixed(a, beta1)
      cm1<-(s1/(2*a)+vy)*(1-exp(-2*a*ta))*exp(-a*tij)
      V <- cm1+(s1*ta*cm2) + na.exclude(me.response) + obs_var_con - mcov - mcov.fixed
    }
  }
  
  ## Function for regression & grid search
  slouch.regression <- function(hl_vy){
    hl <- hl_vy[1]; vy <- hl_vy[2]
    
    if(hl == 0){
      a <- Inf
    }else{
      a <- log(2)/hl
      cm2 <- make.cm2(a,tia,tja,ta,N,T.term)
    }
    X <- make.X(a)
    
    ## beta1 exists from OLS estimate
    
    con.count <- 0
    repeat{
      V <- calc.V(hl, vy, beta1, cm2)
      
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
    gls.resid<-Y-eY
    
    log.det.V <- mk.log.det.V(V = V, N = N)
    
    sup1 <- -N/2*log(2*pi)-0.5*log.det.V-0.5*(t(resid) %*% V.inverse%*%resid)
    print(as.numeric(round(cbind(hl, vy, sup1, t(beta1)), 4)))
    
    list(support = sup1,
         V = V,
         beta1 = beta1,
         X = X,
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








make.slouch.regression <- function(){
 
}