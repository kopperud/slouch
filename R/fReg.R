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
  
  obs_var_con <- mk.obs_var_con(a, hl, beta1, T, N, xx, x.ols, error_condition) ## Putting obs_var_con outside iterated loop works. Why?
  
  ## Initial "half"-V
  V.initial <- ((vy)*(1-exp(-2*a*ta))*exp(-a*tij)) + na.exclude(me.response)
  
  ##### iterated GLS
  con.count<-0;  # Counter for loop break if Beta's dont converge #
  repeat
  {
    #obs_var_con <- mk.obs_var_con(a, hl, beta1, T, N, xx, x.ols, error_condition) ## This causes big trouble for beta1 to converge, only when hl == 0. No idea why.
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











make.seed.fReg.ffANCOVA <- function(treepar, modelpar){
  list2env(treepar, envir = environment())
  list2env(modelpar, envir = environment())
  
  if(model.type=="fReg"){
    x.ols<-cbind(1, fixed.pred)
    beta1<-solve(t(x.ols)%*%x.ols)%*%(t(x.ols)%*%Y)
    n.fixed<-1
    regime.specs <- NULL
  }
  else
    {
    
    regime.specs<-fixed.fact;
    n.fixed<-length(levels(as.factor(regime.specs)))
    regime.specs<-as.factor(regime.specs)

    #treepar$regime.specs <- regime.specs
    #print(regime.specs)
    x.ols<-weight.matrix(10, topology, times, N, regime.specs, fixed.pred, intercept);


    for (i in length(x.ols[1,]):1){
      if(sum(x.ols[,i]) == 0) {
        x.ols<-x.ols[,-i]; n.fixed<-n.fixed-1
        }  #removes internal regimes that have only zero entries
    }
    beta1<-solve(t(x.ols)%*%x.ols)%*%(t(x.ols)%*%Y)
  }
  
  
  if(model.type=="fReg"){
    Vu<-diag(c(rep(0,N), c(as.numeric(na.exclude(me.fixed.pred)))))
  }else{
    Vu<-diag(c(rep(0,n.fixed*N),as.numeric(na.exclude(me.fixed.pred))))
  }
  
  #print(Vu)
  
  true_var<-matrix(data=0, ncol=n.fixed.pred, nrow=N)
  for (i in 1:n.fixed.pred){
    true_var[,i]<-var(na.exclude(fixed.pred[,i]))-as.numeric(na.exclude(me.fixed.pred[,i]))
  }
  true_var<-c(true_var)
  
  
  if(model.type=="fReg") {
    Vd<-diag(c(rep(0,N),true_var))
  }  else {
    Vd<-diag(c(rep(0,n.fixed*N), true_var))
  }
  
  error_condition <- Vu - (Vu%*%pseudoinverse(Vu+Vd)%*%Vu)
  
  ## Multiplies with betas ##
  xx<-seq(from=1, to=length(Vu[,1]), by=N)
  
  
  if(length(xx) > length(beta1)) {
    beta1<-as.matrix(c(0, beta1))
  }
  
  
  obs_var_con <-matrix(0, nrow=N, ncol=N)
  for (e in seq(from=1, to=ncol(x.ols), by=1)){
    for (j in seq(from=1, to=ncol(x.ols), by=1)) {
      tmp<-error_condition[xx[e]:(e*N),xx[j]:(j*N)]*beta1[e]*beta1[j]
      obs_var_con <-obs_var_con + tmp
    }
  }
  
  #beta1<-solve(t(x.ols)%*%x.ols)%*%(t(x.ols)%*%Y)
  
  list(x.ols = x.ols,
       xx = xx,
       Vd = Vd,
       Vu = Vu,
       beta1 = beta1,
       regime.specs = regime.specs,
       error_condition = error_condition,
       obs_var_con = obs_var_con,
       n.fixed = n.fixed)
  
}