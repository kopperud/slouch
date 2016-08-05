IntcptReg <- function(hl_vy, treepar,modelpar){
  list2env(treepar, envir = environment())
  list2env(modelpar, envir = environment())
  
  hl <- hl_vy[1]; vy <- hl_vy[2]
  if(hl==0)
  {
    a<-1000000000000000000000;
    V<-diag(rep(vy, times=N)) + me.response;
  }
  else
  {
    a <- log(2)/hl;
    V<-((vy)*(1-exp(-2*a*ta))*exp(-a*tij))+me.response;  # The V matrix will get singular when vy and h are 0, when the model does not contain any measurement error.
  }
  if(model.type=="IntcptReg")
  {
    if(hl==0 || a>=1000000000000000000000) {
      X<-matrix(data=1, nrow=N, ncol=1)
    }
    else
      if(ultrametric==TRUE) {
        X<-matrix(data=1, nrow=N, ncol=1)
      } 
      else
      {
        X<-matrix(data=0, nrow=N, ncol=2)
        X[,1]<-1-exp(-a*T)
        X[,2]<-exp(-a*T)
      }
  }
  else{
    X<-weight.matrix(a, topology,times, N, regime.specs, fixed.cov, intercept)
  }

  # GLS estimation of parameters for fixed model
  V.inverse<-solve(V)

  tmp<-pseudoinverse(t(X)%*%V.inverse%*%X) #### Ask Thomas about this one. Bjorn: "tmp" means variance of beta.i ? "beta.i.var ????
  if(Inf %in% tmp) {print("Pseudoinverse of (XT V?X)?1 contained values = Inf, which were set to 10^300")}
  tmp <-replace(tmp, tmp ==Inf, 10^300)
  if(-Inf %in% tmp) {print("Pseudoinverse of (XT V?X)?1 contained values = -Inf, which were set to -10^300")}
  tmp <-replace(tmp, tmp ==-Inf, -10^300)

  beta.i<-tmp%*%(t(X)%*%V.inverse%*%Y)
  beta0<-beta.i
  eY<-X%*%beta0
  resid<-Y-eY
  log.det.V <- mk.log.det.V(V = V, N = N)
  
  sup1 <- -N/2*log(2*pi)-0.5*log.det.V-0.5*(t(resid) %*% V.inverse%*%resid)
  print(as.numeric(round(cbind(if(a!=0)log(2)/a else 0.00, vy, sup1, t(beta0)), 4)))
  if(sup1 == "NaN") sup1 <- -10^100
  #return(sup1)
  
  ## Return list of variables
  list(support = sup1,
       V = V,
       beta0 = beta0,
       X = X,
       beta.i.var = tmp,
       alpha.est = a,
       vy.est = vy)
}
