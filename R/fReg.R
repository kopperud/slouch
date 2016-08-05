## Function to fit ffANCOVA + fReg

fReg <- function(hl_vy, treepar, modelpar, seed){
  list2env(modelpar, envir = environment())
  list2env(treepar, envir = environment())
  list2env(seed, envir = environment())
  
  hl <- hl_vy[1]; vy <- hl_vy[2]
  if(hl==0)
  {
    a<-1000000000000000000000
  }
  else
  {
    a <- log(2)/hl
  }
  
  if(model.type=="fReg"){
    X<-cbind(1, fixed.pred)
  }else{
    X<-weight.matrix(a, topology, times, N, regime.specs, fixed.pred, intercept)
  }
  
  ## Initial "half"-V
  V.initial <- ((vy)*(1-exp(-2*a*ta))*exp(-a*tij)) + na.exclude(me.response)
  
  ##### iterated GLS
  con.count<-0;  # Counter for loop break if Beta's dont converge #
  repeat
  {
    mk.obs_var_con(a, hl, beta1, T, N, xx, x.ols, error_condition)
    V <- V.initial + obs_var_con - diag(as.vector(na.exclude(me.cov%*%(2*beta1[(n.fixed+1):length(beta1),]))))

    # INTERMEDIATE ESTIMATION OF OPTIMAL REGRESSION #
    V.inverse<-solve(V)
    beta.i.var <- pseudoinverse(t(X)%*%V.inverse%*%X)
    beta.i<-beta.i.var%*%(t(X)%*%V.inverse%*%Y)
    
    if(length(beta.i) > length(beta1)) { ##### What is this ???????? Why on earth would beta.i & beta1 not have same dimensions ??? Bjorn
      beta1<-as.matrix(c(0, beta1))
    }
    
    con.count <- con.count + 1
    if (test.conv.fReg(beta.i = beta.i, beta1 = beta1, convergence = convergence, con.count = con.count, ultrametric = ultrametric)) {
      break
    }
    beta1<-beta.i
  }
  
  eY<-X%*%beta1
  resid<-Y-eY
  
  log.det.V <- mk.log.det.V(V = V, N = N)
  
  #gof[i, k] <- -N/2*log(2*pi)-0.5*log.det.V-0.5*(t(resid) %*% V.inverse%*%resid);
  sup1 <- -N/2*log(2*pi)-0.5*log.det.V-0.5*(t(resid) %*% V.inverse%*%resid)
  print(as.numeric(round(cbind(if(a!=0)log(2)/a else 0.00, vy, sup1, t(beta1)), 4)))
  
  ## Return list of variables
  list(support = sup1,
       V = V,
       beta1 = beta1,
       X = X,
       beta.i.var = beta.i.var,
       alpha.est = a,
       vy.est = vy)
}