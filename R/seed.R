## Function to seed the OLS
## outputs; beta1, x.ols, xx, yy, Vu, Vd
ols.seed <- function(treepar, modelpar){
  list2env(treepar, envir = environment())
  list2env(modelpar, envir = environment())
  
  ## RANDOM PREDICTOR THETA AND SIGMA ESTIMATES
  
  
  
  if(!is.null(random.cov)){
    n.pred<-length(as.matrix(random.cov)[1,])
    pred<-data.frame(random.cov)
    pred<-matrix(data=pred[!is.na(pred)], ncol=n.pred)
    
    if(!is.null(me.random.cov)){
      me.pred<-matrix(data=me.random.cov[!is.na(me.random.cov)], ncol=n.pred)
    }else{
      me.pred<-matrix(data=0, nrow=N, ncol=n.pred)
    }
    
    s.X<-matrix(data=0, ncol=n.pred)  # PREDICTOR SIGMA
    theta.X<-matrix(data=0, ncol=n.pred)  #PREDICTOR THETA
    for(i in 1:n.pred)
    {
      s.X[,i] <- as.numeric(sigma.X.estimate(pred[,i],me.pred[,i], topology, times)[2])
      theta.X[,i] <- as.numeric(sigma.X.estimate(pred[,i],me.pred[,i], topology, times)[1])
    }
  }else{
    n.pred <- 0
    pred <- NULL
    s.X <- NULL
    theta.X <- NULL
    me.pred<-matrix(data=0, nrow=N, ncol=n.pred)
  }
  
  if(is.null(mecov.random.cov)){
    me.cov<-matrix(data=0, nrow=N, ncol=n.pred)
  }else{
    me.cov<-matrix(data=mecov.random.cov[!is.na(mecov.random.cov)], ncol=ncol(as.matrix(mecov.random.cov)))
  }
  # END OF RANDOM PREDICTOR THETA AND SIGMA ESTIMATES
  
  # FIXED COVARIATES
  if(!is.null(fixed.cov)){
    n.fixed.pred<-length(as.matrix(fixed.cov)[1,])
    
    fixed.pred<-data.frame(fixed.cov)
    fixed.pred<-matrix(data=fixed.pred[!is.na(fixed.pred)], ncol=n.fixed.pred)
    
    
  }else{
    n.fixed.pred <- 0
    fixed.pred <- NULL
  }
  
  if(is.null(me.fixed.cov)) {
    me.fixed.pred<-matrix(data=0, nrow=N, ncol=n.fixed.pred)
  } else{
    me.fixed.pred<- matrix(data=me.fixed.cov[!is.na(me.fixed.cov)], ncol=n.fixed.pred)
  }
  
  if(is.null(mecov.fixed.cov)){
    me.fixed.cov<-matrix(data=0, nrow=N, ncol=n.fixed.pred)
  }else{
    me.fixed.cov<-matrix(data=mecov.fixed.cov[!is.na(mecov.fixed.cov)], ncol=n.fixed.pred)
  }
  
  regime.specs<-as.factor(fixed.fact)
  n.factor<-length(levels(regime.specs))
  
  if(!is.null(fixed.fact) == TRUE){
    #x.ols<-cbind(weight.matrix(10, topology, times, N, regime.specs, fixed.pred, intercept), pred)
    x.ols<-cbind(weight.matrix(10, topology, times, N, regime.specs, fixed.cov, intercept), pred)
  }else{
    x.ols<-matrix(cbind(1, fixed.pred, pred), nrow=N)
  }
  beta1<-solve(t(x.ols)%*%x.ols)%*%(t(x.ols)%*%Y)
  #if(ultrametric == FALSE & !is.null(random.cov)){
  if(is.null(intercept)){
    beta1<-rbind(0, 0, beta1) # 2 additional parameter seeds for Ya and Xa
  }
  
  
  #Defining the dimensionality of Vd ## Might need to change the dimensionality of this matrix for different versions of the model and how that impacts the number of columns in x.ols
  Vd<-matrix(0,ncol=(N*length(beta1[,1])), nrow=(N*length(beta1[,1])))
  
  
  
  
  
  
  ## Define true_var
  if(!is.null(fixed.pred)){
    #Putting in elements in VD for fixed covariates
    true_var<-matrix(data=0, ncol=n.fixed.pred, nrow=N);
    for (i in 1:n.fixed.pred)
    {
      true_var[,i]<-var(na.exclude(fixed.pred[,i]))-as.numeric(na.exclude(me.fixed.pred[,i]))
    }
    true_var<-c(true_var)
    
    if(is.null(random.cov)){
      if(is.null(fixed.fact)){
        Vd<-diag(c(rep(0,N),true_var))
      }else{
        Vd<-diag(c(rep(0,n.factor*N), true_var))
      }
    }else{
      Vd[(((N* n.factor)+(1)):(((N* n.factor)+(1))+((N*n.fixed.pred)-1))),(((N* n.factor)+(1)):(((N* n.factor)+(1))+((N*n.fixed.pred)-1)))]<-diag(c(true_var))
    }
    
    # if(model.type=="fReg") {
    #   Vd<-diag(c(rep(0,N),true_var))
    # }
    # if(model.type =="ffANCOVA"){
    #   Vd<-diag(c(rep(0,n.factor*N), true_var))
    # }
    # 
    # if(model.type =="mmfANCOVA" | model.type =="mfReg"){
    #   Vd[(((N* n.factor)+(1)):(((N* n.factor)+(1))+((N*n.fixed.pred)-1))),(((N* n.factor)+(1)):(((N* n.factor)+(1))+((N*n.fixed.pred)-1)))]<-diag(c(true_var))
    # }
  }
  
  
  #Putting in elements in VD for random covariates
  xx<-seq(from=1, to = length(Vd[,1]), by=N)
  yy<-seq(from=N, to = length(Vd[,1]), by=N)
  
  if(ultrametric == TRUE){
    xx<-xx[-(1:(n.factor+n.fixed.pred))]
    yy<-yy[-(1:(n.factor+n.fixed.pred))]
  }else{
    xx<-xx[-(1:(2 + n.factor+n.fixed.pred))]
    yy<-yy[-(1:(2 + n.factor+n.fixed.pred))]
  }
  
  if(!is.null(random.cov)){
    for (i in seq(from=1, to=nrow(s.X), by=1)){
      Vd[xx[i]:yy[i],xx[i]:yy[i]] <- pt$bt*s.X[,i]
    }
  }
  
  ## ------------------------------------------ ##
  ##               Defining Vu                  ##
  ## ------------------------------------------ ##
  
  if(!is.null(random.cov) & !is.null(fixed.fact)){
    if(ultrametric == TRUE){
      Vu<-diag(c(rep(0,(N*(n.factor))), 
                 c(as.numeric(na.exclude(me.fixed.pred))),
                 c(as.numeric(na.exclude(me.pred)))))
    }  else{
      Vu<-diag(c(rep(0,N*(2+ n.factor))), 
               c(as.numeric(na.exclude(me.fixed.pred))), 
               c(as.numeric(na.exclude(me.pred))))
    } 
  }
  if(!is.null(random.cov) & is.null(fixed.cov) & is.null(fixed.fact)){
    if(ultrametric == TRUE) {
      Vu<-diag(c(rep(0,N), c(as.numeric(na.exclude(me.pred)))))
    } else{
      Vu<-diag(c(rep(0,N*3), c(as.numeric(na.exclude(me.pred)))))
    } 
  }
  
  if(is.null(random.cov) & !is.null(fixed.cov) & !is.null(fixed.fact)){
    Vu<-diag(c(rep(0,n.factor*N), as.numeric(na.exclude(me.fixed.pred))))
  }
  
  if(is.null(random.cov) & !is.null(fixed.cov) & is.null(fixed.fact)){
    Vu<-diag(c(rep(0,N), c(as.numeric(na.exclude(me.fixed.pred)))))
  }
  
  if(is.null(random.cov) & is.null(fixed.cov)){
    Vu <- NULL
  }
  
  if(!is.null(random.cov) | !is.null(fixed.cov)){
    xx<-seq(from=1, to=length(Vu[,1]), by=N)
    error_condition<-Vu-(Vu%*%pseudoinverse(Vu+Vd)%*%Vu)
  }else{
    xx <- NULL
    error_condition <- NULL
  }

  
  
  
  list(n.pred = n.pred,
       s.X = s.X,
       theta.X = theta.X,
       pred = pred,
       me.pred = me.pred,
       me.cov = me.cov,
       n.factor = n.factor,
       fixed.pred = fixed.pred,
       n.fixed.pred = n.fixed.pred,
       me.fixed.pred = me.fixed.pred,
       me.fixed.cov = me.fixed.cov,
       x.ols = x.ols,
       ols.beta1 = beta1,
       Y = Y,
       Vd = Vd,
       Vu = Vu,
       xx = xx,
       yy = yy,
       error_condition = error_condition)
}