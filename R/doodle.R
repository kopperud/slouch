regression.closures <- function(treepar, modelpar, seed){
  ## bind global variables to this environment
  list2env(treepar, envir = environment())
  list2env(modelpar, envir = environment())
  list2env(seed, envir = environment())
  

  
  ## Establish function to calculate design matrix X
  if(is.null(fixed.fact) == FALSE){
    calc.X <- function(a, hl){
      if(model.type == "ffANOVA"){
        matrix(cbind(1-exp(-a*T.term), exp(-a*T.term)), nrow=N)
      }else{
        
      }
      if(hl == 0){
        cbind(weight.matrix(a, topology, times, N, regime.specs, fixed.pred, intercept), pred)
      }else{
        cbind(weight.matrix(a, topology, times, N, regime.specs, fixed.pred, intercept), (1-(1-exp(-a*T.term))/(a*T.term))*pred)
      }
    }
  }else{
    calc.X <- function(a, hl){
      if(hl == 0){
        cbind(1, fixed.pred, pred)
      }else{
        if (ultrametric == TRUE){
          matrix(cbind(1, fixed.pred, (1-(1-exp(-a*T.term))/(a*T.term))*pred), nrow=N)
        }else{
          cbind(1-exp(-a*T.term), 1-exp(-a*T.term)-(1-(1-exp(-a*T.term))/(a*T.term)), exp(-a*T.term), fixed.pred, (1-(1-exp(-a*T.term))/(a*T.term))*pred)
        }
      }
    }
  }
  
  
  ## Function to calculate variance of instantaneous predictor in the OU-process
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
    if(!is.null(me.fixed.cov)){
      if(ultrametric == TRUE){
        diag(as.numeric(me.fixed.cov%*%(2*beta1[2:(length(beta1)-n.pred),])))
      }else{
        diag(as.numeric(me.fixed.cov%*%(2*beta1[4:(length(beta1)-n.pred),])))
      }
    }else{
      NULL
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
      s1 <- calc.s1(beta1)
      mcov <- calc.mcov(a, beta1)
      mcov.fixed <- calc.mcov.fixed(a, beta1)
      cm1<-(s1/(2*a)+vy)*(1-exp(-2*a*ta))*exp(-a*tij)
      # for (i in c(cm1, (s1*ta*cm2), na.exclude(me.response), obs_var_con, mcov, mcov.fixed)){
      #   print(str(i))
      # }
      print(str(cm1))
      print(str(s1*ta*cm2))
      print(str(na.exclude(me.response)))
      print(str(obs_var_con))
      print(str(mcov))
      print(str(mcov.fixed))
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
    X <- calc.X(a, hl)
    
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
  #return(all.closures)
}


## Function to seed the OLS
## outputs; beta1, x.ols, xx, yy, Vu, Vd
ols.seed <- function(treepar, modelpar){
  list2env(treepar, envir = environment())
  list2env(modelpar, envir = environment())
  
  ## RANDOM PREDICTOR THETA AND SIGMA ESTIMATES
  n.pred<-length(as.matrix(random.cov)[1,])
  
  if(!is.null(random.cov)){
    pred<-data.frame(random.cov)
    pred<-matrix(data=pred[!is.na(pred)], ncol=n.pred)
    
    if(!is.null(me.random.cov)){
      me.pred<-matrix(data=me.random.cov[!is.na(me.random.cov)], ncol=n.pred)
    }else{
      me.pred<-matrix(data=0, nrow=N, ncol=n.pred)
    }


    
    if(is.null(mecov.random.cov)){
      me.cov<-matrix(data=0, nrow=N, ncol=n.pred)
    }else{
      me.cov<-matrix(data=mecov.random.cov[!is.na(mecov.random.cov)], ncol=n.pred)
    } 
    
    s.X<-matrix(data=0, ncol=n.pred)  # PREDICTOR SIGMA
    theta.X<-matrix(data=0, ncol=n.pred)  #PREDICTOR THETA
    for(i in 1:n.pred)
    {
      s.X[,i] <- as.numeric(sigma.X.estimate(pred[,i],me.pred[,i], topology, times)[2])
      theta.X[,i] <- as.numeric(sigma.X.estimate(pred[,i],me.pred[,i], topology, times)[1])
    }
  }else{
    pred <- NULL
    s.X <- NULL
    theta.X <- NULL
    me.pred<-matrix(data=0, nrow=N, ncol=n.pred)
    me.cov <- NULL
  }

  
  
  # END OF RANDOM PREDICTOR THETA AND SIGMA ESTIMATES
  
  # FIXED COVARIATES
  n.fixed.pred<-length(fixed.cov[1,])
  if(!is.null(fixed.cov)){
    fixed.pred<-data.frame(fixed.cov)
    fixed.pred<-matrix(data=fixed.pred[!is.na(fixed.pred)], ncol=n.fixed.pred)
    
    if(is.null(me.fixed.cov)) {
      me.fixed.pred<-matrix(data=0, nrow=N, ncol=n.fixed.pred)
    } else{
      me.fixed.pred<- matrix(data=me.fixed.cov[!is.na(me.fixed.cov)], ncol=n.fixed.pred)
    }
    if(is.null(me.cov.fixed.cov)){
      me.fixed.cov<-matrix(data=0, nrow=N, ncol=n.fixed.pred)
    }else{
      me.fixed.cov<-matrix(data=me.cov.fixed.cov[!is.na(me.cov.fixed.cov)], ncol=n.fixed.pred);
    }
  }else{
    fixed.pred <- NULL
    me.fixed.pred <- matrix(data=0, nrow=N, ncol=n.fixed.pred)
    me.fixed.cov <- NULL
  }
  
  regime.specs<-as.factor(fixed.fact)
  n.factor<-length(levels(regime.specs))
  
  if(!is.null(fixed.fact) == TRUE){
    x.ols<-cbind(weight.matrix(10, topology, times, N, regime.specs, fixed.pred, intercept), pred)
  }else{
    x.ols<-matrix(cbind(1, fixed.pred, pred), nrow=N)
  }
  
  beta1<-solve(t(x.ols)%*%x.ols)%*%(t(x.ols)%*%Y)
  if(ultrametric == FALSE){
    beta1<-rbind(0, 0, beta1) # 2 additional parameter seeds for Ya and Xa
  }
  
  
  #Defining the dimensionality of Vd ## Might need to change the dimensionality of this matrix for different versions of the model and how that impacts the number of columns in x.ols
  Vd<-matrix(0,ncol=(N*length(beta1[,1])), nrow=(N*length(beta1[,1])))
  

  
  if(model.type=="fReg") {
    Vd<-diag(c(rep(0,N),true_var))
  }
  if(model.type =="ffANCOVA"){
    Vd<-diag(c(rep(0,n.factor*N), true_var))
  }
  
  
  ## Define true_var
  if(!is.null(fixed.pred)){
    #Putting in elements in VD for fixed covariates
    true_var<-matrix(data=0, ncol=n.fixed.pred, nrow=N);
    for (i in 1:n.fixed.pred)
    {
      true_var[,i]<-var(na.exclude(fixed.pred[,i]))-as.numeric(na.exclude(me.fixed.pred[,i]))
    }
    true_var<-c(true_var)
    Vd[(((N* n.factor)+(1)):(((N* n.factor)+(1))+((N*n.fixed.pred)-1))),(((N* n.factor)+(1)):(((N* n.factor)+(1))+((N*n.fixed.pred)-1)))]<-diag(c(true_var))
  }

  
  #Putting in elements in VD for random covariates
  
  xx<-seq(from=1, to	=length(Vd[,1]), by=N)
  yy<-seq(from=N, to	=length(Vd[,1]), by=N)
  
  # if(ultrametric == TRUE){
  #   xx<-xx[-(1:(n.factor+n.fixed.pred))]
  #   yy<-yy[-(1:(n.factor+n.fixed.pred))]
  # }  else{
  #   xx<-xx[-(1:(2 + n.factor+n.fixed.pred))]
  #   yy<-yy[-(1:(2 + n.factor+n.fixed.pred))]
  # } 
  
  if(!is.null(random.cov)){
    if(!is.null(random.cov)){
      for (i in seq(from=1, to=nrow(s.X), by=1)){
        Vd[xx[i]:yy[i],xx[i]:yy[i]]<-pt$bt*s.X[,i]
      }
    }
  }

  
  # Defining Vu
  if(!is.null(random.cov) & !is.null(fixed.cov)){
    if(ultrametric == TRUE){
      Vu<-diag(c(rep(0,(N*(n.factor))), c(as.numeric(na.exclude(me.fixed.pred))),c(as.numeric(na.exclude(me.pred)))))
    }  else{
      Vu<-diag(c(rep(0,N*(2+ n.factor))), c(as.numeric(na.exclude(me.fixed.pred))), c(as.numeric(na.exclude(me.pred))))
    } 
  }
  if(!is.null(random.cov) & is.null(fixed.cov)){
    if(ultrametric == TRUE) {
      Vu<-diag(c(rep(0,N), c(as.numeric(na.exclude(me.pred)))))
      print(sum(Vu));print(sum(Vu))
    } else{
      Vu<-diag(c(rep(0,N*3), c(as.numeric(na.exclude(me.pred)))))
    } 
  }
  if(!is.null(fixed.cov) & is.null(random.cov)){
    if(is.null(fixed.fact)){
      Vu<-diag(c(rep(0,N), c(as.numeric(na.exclude(me.fixed.pred)))))
    }else{
      Vu<-diag(c(rep(0,n.fixed*N),as.numeric(na.exclude(me.fixed.pred))))
    }
  }
  
  

  
  error_condition<-Vu-(Vu%*%pseudoinverse(Vu+Vd)%*%Vu)
  

  
  list(n.pred = n.pred,
       s.X = s.X,
       theta.X = theta.X,
       pred = pred,
       me.pred = me.pred,
       me.cov = me.cov,
       fixed.pred = fixed.pred,
       n.fixed.pred = n.fixed.pred,
       me.fixed.pred = me.fixed.pred,
       me.fixed.cov = me.fixed.cov,
       x.ols = x.ols,
       beta1 = beta1,
       Vd = Vd,
       Vu = Vu,
       xx = xx,
       yy = yy,
       error_condition = error_condition)
}
