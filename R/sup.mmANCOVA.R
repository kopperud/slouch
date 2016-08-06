sup.mmANCOVA <- function(hl_vy, N, me.response, ta, tij, T, topology, times, model.type, ultrametric, Y, fixed.cov, pred, xx, beta1, error_condition, s.X, n.pred, num.prob, tia, tja, cm2, me.pred, me.cov, convergence, n.fixed, x.ols, regime.specs, intercept){
  hl <- hl_vy[1]; vy <- hl_vy[2]

  if(hl==0)
  {
    a<-1000000000000000000000
    pred.coef <- 1
  }
  else
  {
    a <- log(2)/hl
    pred.coef <- (1-(1-exp(-a*T))/(a*T))
  }
  X<-cbind(weight.matrix(a, topology, times, N, regime.specs, fixed.cov=NULL, intercept), pred.coef*pred)

  if(length(X[1,]) > length(beta1)) {beta1<-as.matrix(c(0, beta1)); n.fixed<-n.fixed+1 ;}
  if(length(X[1,]) < length(beta1)) {beta1<-solve(t(x.ols)%*%x.ols)%*%(t(x.ols)%*%Y);n.fixed<-length(levels(as.factor(regime.specs))); print("The Ya parameter is dropped as its coefficient is too small");}

  # CODE FOR ESTIMATING BETA USING ITERATED GLS
  con.count<-0;  # Counter for loop break if Beta's dont converge #
  repeat
  {
    obs_var_con <- mk.obs_var_con(a, hl, beta1, T, N, xx, x.ols, error_condition)
    if(hl==0)
    {
      V<-diag(rep(vy, times=N))+na.exclude(me.response)+ obs_var_con - diag(as.numeric(me.cov%*%(2*beta1[(n.fixed+1):length(beta1),])))
    }
    else
    {
      # obs_var_con <- mk.obs_var_con(a, hl, beta1, T, N, xx, x.ols, error_condition)
      s1 <- as.numeric(s.X%*%(beta1[(n.fixed+1):length(beta1),]*beta1[(n.fixed+1):length(beta1),]));
      cm1 <- (s1/(2*a)+vy)*(1-exp(-2*a*ta))*exp(-a*tij)
      cm2 <- make.cm2(a,tia,tja,ta,N,T)
      mcov <- diag(rowSums(matrix(data=as.numeric(me.cov)*t(kronecker(2*beta1[(n.fixed+1):length(beta1),], (1-(1-exp(-a*T))/(a*T)))), ncol=n.pred)));
      V<-cm1 + (s1*ta*cm2) + na.exclude(me.response) + obs_var_con - mcov
    } # END OF If ELSE CONDITION FOR HALF-LIFE 0 OR NOT
    # INTERMEDIATE ESTIMATION OF OPTIMAL REGRESSION #

    V.inverse<-solve(V)
    if(hl==0)
    {
      for (z in length(X[1,]):1)
      {
        if(sum(X[,z]) == 0) {X<-X[,-z]; n.fixed-z}  #removes internal regimes that have only zero entries
      }
    }
    beta.i.var <- pseudoinverse(t(X)%*%V.inverse%*%X)
    beta.i<-beta.i.var%*%(t(X)%*%V.inverse%*%Y)

    con.count <- con.count + 1
    if (test.conv.rReg(beta.i = beta.i, beta1 = beta1, n.pred = n.pred, convergence = convergence, con.count = con.count, ultrametric = ultrametric)) {
      break
    }
    beta1<-beta.i
  }

  V.inverse<-solve(V)
  eY<-X%*%beta1
  resid<-Y-eY

  log.det.V <- mk.log.det.V(V = V, N = N)
  sup1 <- -N/2*log(2*pi)-0.5*log.det.V-0.5*(t(resid) %*% V.inverse%*%resid)
  print(as.numeric(round(cbind(if(a!=0)log(2)/a else 0.00, vy, sup1, t(beta1)), 4)))

  list(support = sup1,
       V = V,
       beta1 = beta1,
       X = X,
       beta1.var = beta.i.var,
       alpha.est = a,
       vy.est = vy)
}
