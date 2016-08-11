## Function to regress the mfReg

sup.mfReg <- function(hl_vy, N, me.response, ta, tij, T.term, topology, times, model.type, ultrametric, Y, fixed.cov, pred, xx, beta1, error_condition, s.X, n.pred, num.prob, tia, tja, cm2, me.pred, me.cov, convergence, n.fixed, fixed.pred, n.fixed.pred, obs_var_con, me.fixed.cov, me.fixed.pred){
  ## Split hl_vy vector
  hl <- hl_vy[1]; vy <- hl_vy[2]
  x.ols<-cbind(1,fixed.pred, pred)
  
  ## Initial OLS estimates
  if(hl==0)
  {
    a<-Inf
    beta1<-solve(t(x.ols)%*%x.ols)%*%(t(x.ols)%*%Y)
  }
  else
  {
    a <- log(2)/hl
    T.term <- T.term
    cm2 <- make.cm2(a,tia,tja,ta,N,T.term)
    if(ultrametric==TRUE)
      beta1<-solve(t(x.ols)%*%x.ols)%*%(t(x.ols)%*%Y)
    else
      beta1<-rbind(0, 0, solve(t(x.ols)%*%x.ols)%*%(t(x.ols)%*%Y))
  }
  
  ### CODE FOR ESTIMATING BETA USING ITERATED GLS ###
  con.count<-0;  # Counter for loop break if Beta's dont converge #
  repeat
  {
    if(hl==0)
    {
      #s1<-as.numeric(s.X%*%(beta1[(2+n.fixed.pred):(n.pred+1+n.fixed.pred),]*beta1[(2+n.fixed.pred):(n.pred+1+n.fixed.pred),]))
      X<-cbind(1, fixed.pred, pred)
      V<-diag(rep(vy, times=N))+na.exclude(me.response)+ obs_var_con - diag(as.numeric(me.cov%*%(2*beta1[(2+n.fixed.pred):(n.pred+n.fixed.pred+1),]))) -diag(as.numeric(me.fixed.cov%*%(2*beta1[2:(length(beta1)-n.pred),])))
    }
    else
    {
      if(ultrametric==TRUE){
        s1<-as.numeric(s.X%*%(beta1[(2+n.fixed.pred):(n.pred+n.fixed.pred+1),]*beta1[(2+n.fixed.pred):(n.pred+n.fixed.pred+1),]))
        X<-cbind(1, fixed.pred, (1-(1-exp(-a*T.term))/(a*T.term))*pred)
        mcov<-diag(rowSums(matrix(data=as.numeric(me.cov)*t(kronecker(2*beta1[(2+n.fixed.pred):(n.pred+n.fixed.pred+1),], (1-(1-exp(-a*T.term))/(a*T.term)))), ncol=n.pred)))
        
        ## fixed.cov measurement error covariances, to be subtracted in V
        last_term <- diag(as.numeric(me.fixed.cov%*%(2*beta1[2:(length(beta1)-n.pred),])))
      }
      else{
        s1<-as.numeric(s.X%*%(beta1[(4+n.fixed.pred):(n.pred+n.fixed.pred+3),]*beta1[(4+n.fixed.pred):(n.pred+n.fixed.pred+3),]))
        X<-cbind(1-exp(-a*T.term), 1-exp(-a*T.term)-(1-(1-exp(-a*T.term))/(a*T.term)), exp(-a*T.term), fixed.pred, (1-(1-exp(-a*T.term))/(a*T.term))*pred)
        mcov<-diag(rowSums(matrix(data=as.numeric(me.cov)*t(kronecker(2*beta1[(4+n.fixed.pred):(n.pred+n.fixed.pred+3),], (1-(1-exp(-a*T.term))/(a*T.term)))), ncol=n.pred)))
        
        ## fixed.cov measurement error as covariances, to be subtracted in V
        last_term <- diag(as.numeric(me.fixed.cov%*%(2*beta1[4:(length(beta1)-n.pred),])))
      }
      cm1<-(s1/(2*a)+vy)*(1-exp(-2*a*ta))*exp(-a*tij)
      
      
      
      V<-cm1+(s1*ta*cm2) + na.exclude(me.response) + obs_var_con - mcov - last_term
    } # END OF ELSE CONDITION FOR HALF-LIFE = 0
    
    # INTERMEDIATE ESTIMATION OF OPTIMAL REGRESSION #
    V.inverse<-solve(V)
    beta.i.var <- pseudoinverse(t(X)%*%V.inverse%*%X)
    beta.i<-beta.i.var%*%(t(X)%*%V.inverse%*%Y)
    
    ## Check for convergence
    if (test.conv(beta.i, beta1, convergence, con.count, ultrametric)) break
    con.cont <- con.count +1 
    
    beta1<-beta.i
  }                            # END OF ITERATED GLS REPEAT LOOP #
  
  ## Use the last beta.i
  beta1<-beta.i
  
  ### END OF ITERATED GLS ESTIMATION FOR BETA #
  
  ## Compute residuals
  eY<-X%*%beta1
  resid<-Y-eY
  
  ## Compute support
  sup1 <- -N/2*log(2*pi)-0.5*log(det(V))-0.5*(t(resid) %*% V.inverse%*%resid)
  
  ## Print hl, vy, support, beta estimates
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

   