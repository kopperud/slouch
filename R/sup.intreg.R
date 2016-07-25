hl.intreg <- function(hl_vy, N, me.response, ta, tij, T, topology, times, regime.specs, model.type, ultrametric, Y, fixed.cov){
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
    if(hl==0 ||a>=1000000000000000000000) X<-matrix(data=1, nrow=N, ncol=1)
    else
      if(ultrametric==TRUE) X<-matrix(data=1, nrow=N, ncol=1)
      else
      {
        X<-matrix(data=0, nrow=N, ncol=2);
        X[,1]<-1-exp(-a*T);
        X[,2]<-exp(-a*T)
      }
  }
  else
    X<-weight.matrix(a, topology,times, N, regime.specs, fixed.cov, intercept)
  # GLS estimation of parameters for fixed model
  V.inverse<-solve(V)

  tmp<-pseudoinverse(t(X)%*%V.inverse%*%X) #### Ask Thomas about this one
  if(Inf %in% tmp) {print("Pseudoinverse of (XT V?X)?1 contained values = Inf, which were set to 10^300")};
  tmp <-replace(tmp, tmp ==Inf, 10^300);
  if(-Inf %in% tmp) {print("Pseudoinverse of (XT V?X)?1 contained values = -Inf, which were set to -10^300")} ;
  tmp <-replace(tmp, tmp ==-Inf, -10^300)

  beta.i<-tmp%*%(t(X)%*%V.inverse%*%Y)
  beta0<-beta.i
  eY<-X%*%beta0
  resid<-Y-eY
  det.V<-det(V)
  if(det.V==0){
    print(paste("Warning: Determinant of V = 0"))
    #Minimum value of diagonal scaling factor
    inv.min.diag.V<-1/min(diag(V))
    V<-V*inv.min.diag.V
    #Rescale and log determinant
    log.det.V<-log(det(V))+log(min(diag(V)))*N
  }
  else {log.det.V<-log(det.V)}
  sup1 <- -N/2*log(2*pi)-0.5*log.det.V-0.5*(t(resid) %*% V.inverse%*%resid);
  if(sup1 == "NaN") sup1 <- -10^100
  return(sup1)
}
