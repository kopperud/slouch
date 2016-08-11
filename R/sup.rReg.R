### Function to return support values for each hl and vy for rReg
# sup.rReg <- function(hl_vy, N, me.response, ta, tij, T.term, topology, times, model.type, ultrametric, Y, fixed.cov, pred, xx, beta1, error_condition, s.X, n.pred, num.prob, tia, tja, cm2, me.pred, me.cov, convergence, n.fixed,make.cm2) {
make.sup.rReg <- function(modelpar,treepar,seed){
  list2env(modelpar, envir = environment())
  list2env(treepar, envir = environment())
  list2env(seed, envir = environment())
  
  function(hl_vy) {  
    
    
    hl <- hl_vy[1]; vy <- hl_vy[2]
    x.ols<-cbind(1, pred)
    
    ## Find starting beta with OLS estimation to seed the iterative GLS
    if (hl != 0 & ultrametric == FALSE){
      beta1 <- rbind(0, 0, solve(t(x.ols)%*%x.ols)%*%(t(x.ols)%*%Y))
    } else{
      beta1 <- solve(t(x.ols)%*%x.ols)%*%(t(x.ols)%*%Y)
    }
    
    ## Set up design matrix X
    if(hl == 0)
    {
      a <- Inf
      X <- cbind(1,pred)
    }
    else
    {
      a <- log(2)/hl
      cm2 <- make.cm2(a,tia,tja,ta,N,T.term)
      cm1.half <- (1-exp(-2*a*ta))*exp(-a*tij)
      if (ultrametric == TRUE)
        X <- cbind(1, (1-(1-exp(-a*T.term))/(a*T.term))*pred)
      else
        X <- cbind(1-exp(-a*T.term), 1-exp(-a*T.term)-(1-(1-exp(-a*T.term))/(a*T.term)), exp(-a*T.term), (1-(1-exp(-a*T.term))/(a*T.term))*pred)
    }
    
    ### CODE FOR ESTIMATING BETA USING ITERATED GLS ###
    con.count<-0;  # Counter for loop break if Beta's dont converge #
    repeat
    {
      V <- estimate.V.rReg(hl, vy, a, ta, tij, T.term, N, xx, x.ols, error_condition, me.response, me.cov, beta1, n.fixed, n.pred, ultrametric, s.X, cm2, cm1.half)
      
      # INTERMEDIATE ESTIMATION OF OPTIMAL REGRESSION #
      V.inverse<-solve(V)
      beta.i.var <- pseudoinverse(t(X)%*%V.inverse%*%X)
      beta.i<-beta.i.var%*%(t(X)%*%V.inverse%*%Y)
      
      con.count <- con.count + 1
      if (test.conv.rReg(beta.i = beta.i, beta1 = beta1, convergence = convergence, n.pred = n.pred, con.count = con.count, ultrametric = ultrametric)) {
        break
      }
      beta1<-beta.i
      # print(beta1)
      #Sys.sleep(0.3)
    }
    
    resid <- Y - (X%*%beta1)
    log.det.V <- mk.log.det.V(V = V, N = N)
    
    sup1 <- -N/2*log(2*pi)-0.5*log.det.V-0.5*(t(resid) %*% V.inverse%*%resid)
    print(as.numeric(round(cbind(hl, vy, sup1, t(beta1)), 4))) # Will increasing the number of digits avoid problems when plotting grid?
    # return(sup1)
    
    list(support = sup1,
         V = V,
         beta1 = beta1,
         X = X,
         beta1.var = beta.i.var,
         alpha.est = a,
         vy.est = vy)
  }
}





estimate.V.rReg <- function(hl, vy, a, ta, tij, T.term, N, xx, x.ols, error_condition, me.response, me.cov, beta1, n.fixed, n.pred, ultrametric, s.X, cm2, cm1.half){
  obs_var_con <- mk.obs_var_con(a, hl, beta1, T.term, N, xx, x.ols, error_condition)

  if (ultrametric == TRUE){
    mcov <- diag(rowSums(matrix(data=as.numeric(me.cov)*t(kronecker(2*beta1[2:(n.pred+1),], (1-(1-exp(-a*T.term))/(a*T.term)))), ncol=n.pred)))
    s1 <- as.numeric(s.X%*%(beta1[2:(n.pred+1),]*beta1[2:(n.pred+1),]))
  }
  else{
    mcov <- diag(rowSums(matrix(data=as.numeric(me.cov)*t(kronecker(2*beta1[4:(n.pred+3),], (1-(1-exp(-a*T.term))/(a*T.term)))), ncol=n.pred)))
    s1 <- as.numeric(s.X%*%(beta1[4:(n.pred+3),]*beta1[4:(n.pred+3),]))
  }

  if(hl == 0)
  {
    diag(rep(vy, times=N)) + na.exclude(me.response) + obs_var_con - diag(as.vector(na.exclude(me.cov%*%(2*beta1[(n.fixed+1):length(beta1),]))))
  }
  else
  {
    cm1<-(s1/(2*a)+vy)*cm1.half

    return(cm1 + (s1*ta*cm2) + na.exclude(me.response) + obs_var_con - mcov)
  } # END OF ELSE CONDITION FOR HALF-LIFE = 0
}


## Seed function for rReg. Initial conditions to feed the proper GLS estimation
seed.rReg <- function(treepar, modelpar){
  list2env(treepar, envir = environment())
  list2env(modelpar, envir = environment())
  
  x.ols<-cbind(1, pred)
  beta1<-solve(t(x.ols)%*%x.ols)%*%(t(x.ols)%*%Y)
  if(ultrametric == FALSE) {
    beta1<-rbind(0, 0, beta1); # 2 additional parameter seeds for Ya and Xa
  }
  n.fixed<-1
  
  ## Setting up the Vu and Vd matrices ##
  Vd<-matrix(0,ncol=(N*length(beta1[,1])), nrow=(N*length(beta1[,1])))
  xx<-seq(from=1, to	=length(Vd[,1]), by=N)
  yy<-seq(from=N, to	=length(Vd[,1]), by=N)
  
  if(ultrametric == TRUE){
    xx<-xx[-1]
    yy<-yy[-1]
  }else{
    xx<-xx[-(1:3)]
    yy<-yy[-(1:3)]
  }
  
  for (i in seq(from=1, to=nrow(s.X), by=1)){
    Vd[xx[i]:yy[i], xx[i]:yy[i]]<-pt$bt*s.X[,i]
    
  }
  
  if(ultrametric == TRUE) {
    Vu<-diag(c(rep(0,N), c(as.numeric(na.exclude(me.pred)))))
  } else{
    Vu<-diag(c(rep(0,N*3), c(as.numeric(na.exclude(me.pred)))))
  } 
  
  error_condition<-Vu-(Vu%*%pseudoinverse(Vu+Vd)%*%Vu)
  
  
  xx<-seq(from=1, to=length(Vu[,1]), by=N)
  
  obs_var_con <-matrix(0, nrow=N, ncol=N)
  for (e in seq(from=1, to=ncol(x.ols), by=1)){
    for (j in seq(from=1, to=ncol(x.ols), by=1)) {
      tmp<-error_condition[xx[e]:(e*N),xx[j]:(j*N)]*beta1[e]*beta1[j]
      obs_var_con <-obs_var_con + tmp
    }
  }
  

  
  seed <- list(x.ols = x.ols,
               beta1 = beta1,
               xx = xx,
               yy = yy,
               Vd = Vd,
               Vu = Vu,
               error_condition = error_condition,
               obs_var_con = obs_var_con,
               n.fixed = n.fixed)
  
  return(seed)
}


