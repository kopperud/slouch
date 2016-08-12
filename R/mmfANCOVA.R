## Function for regression, the mmfANCOVA

reg.mmfANCOVA <- function(hl_vy, N, me.response, ta, tij, T.term, topology, times, model.type, ultrametric, Y, fixed.cov, pred, xx, beta1, error_condition, s.X, n.pred, num.prob, tia, tja, cm2, me.pred, me.cov, convergence, n.fixed, fixed.pred, n.fixed.pred, obs_var_con, me.fixed.cov, me.fixed.pred, regime.specs, intercept, x.ols){
  hl <- hl_vy[1]; vy <- hl_vy[2]
  
  if(hl==0){
    a<-1000000000000000000000
    X<-cbind(weight.matrix(a, topology, times, N, regime.specs, fixed.pred, intercept), pred)
  }else{
    a <- log(2)/hl
    cm2 <- make.cm2(a,tia,tja,ta,N,T.term)
    X<-cbind(weight.matrix(a, topology, times, N, regime.specs, fixed.pred, intercept), (1-(1-exp(-a*T))/(a*T))*pred)
  }
  
  if(length(X[1,]) > length(beta1)) {
    beta1<-as.matrix(c(0, beta1))
    n.fixed<-n.fixed+1
  }
  if(length(X[1,])< length(beta1)) {
    beta1<-solve(t(x.ols)%*%x.ols)%*%(t(x.ols)%*%Y)
    n.fixed<-length(levels(as.factor(regime.specs)))
    print("The Ya parameter is dropped as its coefficient is too small")
  }
  
  # CODE FOR ESTIMATING BETA USING ITERATED GLS
  
  con.count<-0;  # Counter for loop break if Beta's dont converge #
  repeat
  {
    # Update measurement error in X
    obs_var_con <- matrix(0, nrow=N, ncol=N)
    for (e in seq(from=1, to=ncol(x.ols), by=1)){
      for (j in seq(from=1, to=ncol(x.ols), by=1)) {
        tmp <- error_condition[xx[e]:(e*N),xx[j]:(j*N)]*beta1[e]*beta1[j]
        obs_var_con <- obs_var_con + tmp
      }
    }
    
    if(hl==0)
    {
      V<-diag(rep(vy, times=N))+na.exclude(me.response)+ obs_var_con - 
        diag(as.numeric(me.cov%*%(2*beta1[(n.fixed+1+n.fixed.pred):length(beta1),]))) -
        diag(as.numeric(me.fixed.cov%*%(2*beta1[(n.fixed+1):(length(beta1)-n.pred),])))
    }
    else
    {
      s1<-as.numeric(s.X%*%(beta1[(n.fixed+1+n.fixed.pred):length(beta1),]*beta1[(n.fixed+1+n.fixed.pred):length(beta1),]))
      
      cm1<-(s1/(2*a)+vy)*(1-exp(-2*a*ta))*exp(-a*tij)
      
      #print(me.cov)
      ## Calculate measurement error of covariances
      mcov <- diag(rowSums(matrix(data=as.numeric(me.cov)*t(kronecker(2*beta1[(n.fixed+1+n.fixed.pred):length(beta1),], (1-(1-exp(-a*T))/(a*T)))), ncol=n.pred)))
      mcov.fixed <- diag(rowSums(matrix(data=as.numeric(me.fixed.cov)*t(kronecker(2*beta1[(n.fixed+1):(length(beta1)-n.pred),], rep(1, times=N))), ncol=n.fixed.pred)))
      

      
      V<-cm1+(s1*ta*cm2)+na.exclude(me.response) + obs_var_con - mcov - mcov.fixed
    } # END OF If ELSE CONDITION FOR HALF-LIFE 0 OR NOT
    
    # INTERMEDIATE ESTIMATION OF OPTIMAL REGRESSION #
    V.inverse<-solve(V)
    beta.i.var <- pseudoinverse(t(X)%*%V.inverse%*%X)
    beta.i<-beta.i.var%*%(t(X)%*%V.inverse%*%Y)
    
    ## Check for convergence
    if (test.conv(beta.i, beta1, convergence, con.count, ultrametric)) break
    
    beta1 <- beta.i
  }
  ## Remember use the last beta
  beta1 <- beta.i
  
  ### END OF ITERATED GLS ESTIMATION FOR BETA #

  V.inverse<-solve(V)
  eY<-X%*%beta1
  resid<-Y-eY
  
  ## Calculate log-likelihood, print to console
  sup1 <- -N/2*log(2*pi)-0.5*log(det(V))-0.5*(t(resid) %*% V.inverse%*%resid)
  print(as.numeric(round(cbind(hl, vy, sup1, t(beta1)), 4)))
  
  
  ## Return list of outputs
  list(support = sup1,
       V = V,
       beta1 = beta1,
       X = X,
       beta1.var = beta.i.var,
       alpha.est = a,
       vy.est = vy)
}