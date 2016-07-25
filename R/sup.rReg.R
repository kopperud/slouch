### Function to return support values for each hl and vy for rReg
sup.rReg <- function(hl_vy, N, me.response, ta, tij, T, topology, times, model.type, ultrametric, Y, fixed.cov, pred, xx, beta1, error_condition, s.X, n.pred, num.prob, tia, tja, cm2, me.pred, me.cov, convergence, n.fixed,make.cm2) {
  hl <- hl_vy[1]; vy <- hl_vy[2]
  x.ols<-cbind(1, pred)
  beta1 <- make.beta1.rReg(hl, x.ols, Y, ultrametric)

  ## Set up design matrix X
  if(hl==0)
  {
    a <- Inf
    X <- cbind(1,pred)
  }
  else
  {
    a <- log(2)/hl
    cm2 <- make.cm2(a,tia,tja,ta,N,T)
    if (ultrametric == TRUE)
      X <- cbind(1, (1-(1-exp(-a*T))/(a*T))*pred)
    else
      X <- cbind(1-exp(-a*T), 1-exp(-a*T)-(1-(1-exp(-a*T))/(a*T)), exp(-a*T), (1-(1-exp(-a*T))/(a*T))*pred)
  }

  ### CODE FOR ESTIMATING BETA USING ITERATED GLS ###
  con.count<-0;  # Counter for loop break if Beta's dont converge #
  repeat
  {
    V <- estimate.V.rReg(hl, vy, a, ta, tij, T, N, xx, x.ols, error_condition, me.response, me.cov, beta1, n.fixed, n.pred, ultrametric, s.X, cm2)

    # INTERMEDIATE ESTIMATION OF OPTIMAL REGRESSION #
    V.inverse<-solve(V)
    beta.i.var <- pseudoinverse(t(X)%*%V.inverse%*%X)
    beta.i<-beta.i.var%*%(t(X)%*%V.inverse%*%Y)

    con.count <- con.count + 1
    if (test.conv(beta.i = beta.i, beta1 = beta1, convergence = convergence, n.pred = n.pred, con.count = con.count, ultrametric = ultrametric)) {
      break
    }
    beta1<-beta.i
  }

  eY <- X%*%beta1
  resid<-Y-eY
  log.det.V <- mk.log.det.V(V = V, N = N)

  sup1 <- -N/2*log(2*pi)-0.5*log.det.V-0.5*(t(resid) %*% V.inverse%*%resid)
  print(as.numeric(round(cbind(if(a!=0)log(2)/a else 0.00, vy, sup1, t(beta1)), 4))) # Will increasing the number of digits avoid problems when plotting grid?
  # return(sup1)
  list(support = sup1,
       V = V,
       beta1 = beta1,
       X = X,
       beta1.var = beta.i.var)
}

#' @export
make.cm2 <- function(a,tia,tja,ta,N,T){
    T.row <- replicate(N,T)
    T.col <- t(T.row)
    num.prob <- ifelse(ta == 0, 1, (1-exp(-a*ta))/(a*ta))
    return(((1-exp(-a*T.row))/(a*T.row))*((1-exp(-a*T.col))/(a*T.col))-(exp(-a*tia)*(1-exp(-a*T.row))/(a*T.col) + exp(-a*tja)*(1-exp(-a*T.row))/(a*T.row))*(num.prob))
}

test.conv <- function(beta.i = beta.i, beta1 = beta1, convergence = convergence, n.pred = n.pred, con.count = con.count, ultrametric = ultrametric){
  if (ultrametric == TRUE) {
    fstart <- 0
    y <- 1
  }
  else {
    fstart <- 3
    y <- 3
  }
  test<-matrix(nrow=(n.pred+1))
  for(f in (1+fstart):(n.pred+y))
  {
    if(abs(as.numeric(beta.i[f]-beta1[f]))<=convergence) test[(f-fstart)]=0 else test[(f-fstart)]=1
  }
  if(sum(test)==0) return (TRUE)
  if(con.count >= 50)
  {
    message("Warning, estimates did not converge after 50 iterations, last estimates printed out")
    return(TRUE)
  }
  return(FALSE)
}

mk.log.det.V <- function(V, N){
  det.V<-det(V)
  if(det.V==0){
    print(paste("Warning: Determinant of V = 0"))
    #Minimum value of diagonal scaling factor
    inv.min.diag.V<-1/min(diag(V))
    V<-V*inv.min.diag.V
    #Rescale and log determinant
    log(det(V))+log(min(diag(V)))*N
  }
  else {
    log(det.V)
  }
}

mk.obs_var_con <- function(a, hl, beta1, T, N, xx, x.ols, error_condition){
  if (hl == 0) y <- 1 else y <- ((1-(1-exp(-a*T))/(a*T))*(1-(1-exp(-a*T))/(a*T)))
  obs_var_con <- matrix(0, nrow=N, ncol=N)
  for (e in seq(from=1, to=ncol(x.ols), by=1)){
    for (j in seq(from=1, to=ncol(x.ols), by=1)) {
      tmp <- error_condition[xx[e]:(e*N),xx[j]:(j*N)]*beta1[e]*beta1[j]*y
      obs_var_con <- obs_var_con + tmp
    }
  }
  obs_var_con
}

make.beta1.rReg <- function(hl, x.ols, Y, ultrametric){
  if (hl != 0 & ultrametric == FALSE){
    rbind(0, 0, solve(t(x.ols)%*%x.ols)%*%(t(x.ols)%*%Y))
  } else{
    solve(t(x.ols)%*%x.ols)%*%(t(x.ols)%*%Y)
  }
}

estimate.V.rReg <- function(hl, vy, a, ta, tij, T, N, xx, x.ols, error_condition, me.response, me.cov, beta1, n.fixed, n.pred, ultrametric, s.X, cm2){
  obs_var_con <- mk.obs_var_con(a, hl, beta1, T, N, xx, x.ols, error_condition)

  if (ultrametric == TRUE){
    mcov <- diag(rowSums(matrix(data=as.numeric(me.cov)*t(kronecker(2*beta1[2:(n.pred+1),], (1-(1-exp(-a*T))/(a*T)))), ncol=n.pred)))
    s1 <- as.numeric(s.X%*%(beta1[2:(n.pred+1),]*beta1[2:(n.pred+1),]))
  }
  else{
    mcov <- diag(rowSums(matrix(data=as.numeric(me.cov)*t(kronecker(2*beta1[4:(n.pred+3),], (1-(1-exp(-a*T))/(a*T)))), ncol=n.pred)))
    s1 <- as.numeric(s.X%*%(beta1[4:(n.pred+3),]*beta1[4:(n.pred+3),]))
  }

  if(hl==0)
  {
    diag(rep(vy, times=N)) + na.exclude(me.response) + obs_var_con - diag(as.vector(na.exclude(me.cov%*%(2*beta1[(n.fixed+1):length(beta1),]))))
  }
  else
  {
    cm1<-(s1/(2*a)+vy)*(1-exp(-2*a*ta))*exp(-a*tij)
    return(cm1+(s1*ta*cm2)+na.exclude(me.response)+ obs_var_con - mcov)
  } # END OF ELSE CONDITION FOR HALF-LIFE = 0
}
