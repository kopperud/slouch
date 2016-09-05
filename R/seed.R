## ------------------------------------------ ##
##              This is the worst             ##
##         part of the code. 27 aug 2016      ##
## ------------------------------------------ ##



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
    #me.pred <- NULL
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
    me.fixed.pred <- matrix(0, nrow=N, ncol = n.fixed.pred)
  }
  
  if(is.null(me.fixed.cov)) {
    me.fixed.pred<-matrix(data=0, nrow=N, ncol=n.fixed.pred)
  } else{
    me.fixed.pred<- matrix(data=me.fixed.cov[!is.na(me.fixed.cov)], ncol=n.fixed.pred)
  }
  
  if(is.null(mecov.fixed.cov)){
    mecov.fixed.cov<-matrix(data=0, nrow=N, ncol=n.fixed.pred)
  }else{
    mecov.fixed.cov<-matrix(data=mecov.fixed.cov[!is.na(mecov.fixed.cov)], ncol=n.fixed.pred)
  }
  
  regime.specs<-as.factor(fixed.fact)
  n.factor<-length(levels(regime.specs))
  
  if(!is.null(fixed.fact)){
    x.ols<-cbind(weight.matrix(10, topology, times, N, regime.specs, fixed.cov, intercept), pred)
  }else{
    x.ols<-matrix(cbind(1, fixed.pred, pred), nrow=N)
  }
  beta1<-solve(t(x.ols)%*%x.ols)%*%(t(x.ols)%*%Y)
  if(is.null(intercept) & (!is.null(random.cov) | !is.null(fixed.cov))){
    beta1<-rbind(0,
                 0, 
                 beta1) # 2 additional parameter seeds for Ya and Xa ##### Rather, b0 & b1Xa ?
  }
  

  ## ------------------------------------------------------- ##
  ##                                                         ##
  ##                      Defining Vu                        ##
  ##     See equation (10) in Hansen & Bartoszek 2012        ##
  ##                                                         ##
  ## ------------------------------------------------------- ##
  

  if(!is.null(fixed.pred)){
    true_var_matrix <- list()
    for (i in 1:n.fixed.pred)
    {
      true_var_matrix[[i]] <- diag(var(na.exclude(fixed.pred[,i]))-as.numeric(na.exclude(me.fixed.pred[,i])))
    }
  }else{
    true_var_matrix <- NULL
  }
  
  if(!is.null(random.cov)){
    
    for (i in seq(from = 1, to = nrow(s.X), by=1)){
      #true_var_random_cov <- lapply(s.X, function(x) pt$bt*x). Ask thomas. Which is correct? This or the below. Seems to have nearly equivalent outcome, perhaps rounding error.
    }
    true_var_random_cov <- list()
    for (i in 1:n.pred){
      true_var_random_cov[[i]] <- diag(var(na.exclude(pred[,i]))-as.numeric(na.exclude(me.pred[,i])))
      
    }
  }else{
    true_var_random_cov <- NULL
  }
  
  
  ## All variances
  Vu <- c(true_var_matrix, true_var_random_cov)
  
  if(!is.null(random.cov) | !is.null(fixed.cov)){
    Vu_given_x <- list()
    for (i in 1:length(Vu)){
      Vx <- diag(cbind(me.fixed.pred, me.pred)[,i]) + Vu[[i]]
      Vu_given_x[[i]] <- Vu[[i]] - (Vu[[i]] %*% solve(Vx) %*% Vu[[i]])
    }
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
       Vu_given_x = Vu_given_x)
}