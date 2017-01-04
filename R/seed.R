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
      s.X[,i] <- as.numeric(sigma.X.estimate(phy, ta, pred[,i],me.pred[,i])[2]) ## Change var name to squared, variance. Not SE.
      theta.X[,i] <- as.numeric(sigma.X.estimate(phy, ta, pred[,i],me.pred[,i])[1])
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

  ## ------------------------------------------------------- ##
  ##                                                         ##
  ##                Defining Vu, Vd, Vx, Vu|x                ##
  ##     See equation (10) in Hansen & Bartoszek 2012        ##
  ##                                                         ##
  ## ------------------------------------------------------- ##

  if(!is.null(fixed.pred)){
    Vd_fixed <- list()
    Vu_fixed <- list()
    for (i in 1:n.fixed.pred)
    {
      Vu_fixed[[i]] <- diag(na.exclude(me.fixed.pred[,i]))
      Vd_fixed[[i]] <- diag(rep(var(na.exclude(fixed.pred[,i])), N)) - Vu_fixed[[i]]
      if(any(Vd_fixed[[i]] < 0)){
        Vd_fixed[[i]][Vd_fixed[[i]] < 0] <- 0
        warning("Vd contains negative variances, scaled up to 0. Do any of the predictor variables have a larger measurement error than their trait variance?")
      }
    }
  }else{
    Vd_fixed <- NULL
    Vu_fixed <- NULL
  }
  
  if(!is.null(random.cov)){
    Vu_random <- list()
    Vd_random <- lapply(s.X, function(e) ta*e)
    # for (i in seq(from = 1, to = nrow(s.X), by=1)){
    for(i in seq_along(s.X)){
      Vu_random[[i]] <- diag(na.exclude(me.pred[,i]))
    }
  }else{
    Vd_random <- NULL
    Vu_random <- NULL
  }
  
  Vu <- c(Vu_fixed, Vu_random)
  Vd <- c(Vd_fixed, Vd_random)
  
  if(!is.null(random.cov) | !is.null(fixed.cov)){
    Vu_given_x <- list()
    Vx <- list()
    for (i in 1:length(Vd)){
      Vx[[i]] <- Vd[[i]] + Vu[[i]]
      Vu_given_x[[i]] <- Vu[[i]] - (Vu[[i]] %*% solve(Vx[[i]]) %*% Vu[[i]])
    }
  }else{
    Vu_given_x = NULL
    Vx <- NULL
  }
  
  list(s.X = s.X,
       theta.X = theta.X,
       pred = pred,
       me.pred = me.pred,
       n.pred = n.pred,
       me.cov = me.cov,
       fixed.pred = fixed.pred,
       n.fixed.pred = n.fixed.pred,
       me.fixed.pred = me.fixed.pred,
       me.fixed.cov = me.fixed.cov,
       mecov.fixed.cov = mecov.fixed.cov,
       Vu_given_x = Vu_given_x,
       Vu = Vu,
       Vx = Vx,
       Vd = Vd)
}