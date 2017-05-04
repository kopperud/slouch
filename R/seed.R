## ------------------------------------------ ##
##              This is the worst             ##
##         part of the code. 27 aug 2016      ##
## ------------------------------------------ ##


seed <- function(tree, pars, control){
  list2env(tree, envir = environment())
  list2env(pars, envir = environment())
  list2env(control, envir = environment())
  
  n <- length(tree$phy$tip.label)
  
  if(!is.null(random.cov)){
    random.cov <- cbind(random.cov)
    n.pred <- ncol(random.cov)
    
    if(!is.null(me.random.cov)){
      me.random.cov <- cbind(me.random.cov)
    }else{
      me.random.cov <- matrix(data=0, nrow=n, ncol = n.pred)
    }
    brownian <- mapply(function(y, y_me) sigma.X.estimate(tree$phy, ta, y, y_me),
                       y = split(t(random.cov), colnames(random.cov)),
                       y_me = split(t(me.random.cov), colnames(random.cov)),
                       SIMPLIFY = FALSE)
    
    sigma_squared <- sapply(brownian, function(x) x$sigma_squared)
    theta.X <- sapply(brownian, function(x) x$mean)

  }else{
    n.pred <- 0
    random.cov <- NULL
    sigma_squared <- NULL
    theta.X <- NULL
    me.random.cov<-matrix(data=0, nrow=n, ncol=n.pred)
  }
  

  if(!is.null(fixed.cov)){
    n.fixed.cov <- length(as.matrix(fixed.cov)[1,])
    fixed.cov <- as.matrix(fixed.cov, dimnames = list(NULL, names.fixed.cov))
  }else{
    n.fixed.cov <- 0
    fixed.cov <- NULL
    me.fixed.cov <- matrix(0, nrow=n, ncol = n.fixed.cov)
  }
  
  if(is.null(me.fixed.cov)) {
    me.fixed.cov<-matrix(0, nrow=n, ncol=n.fixed.cov)
  } else{
    me.fixed.cov<- matrix(me.fixed.cov, ncol = n.fixed.cov)
  }

  ## ------------------------------------------------------- ##
  ##                                                         ##
  ##                Defining Vu, Vd, Vx, Vu|x                ##
  ##     See equation (10) in Hansen & Bartoszek 2012        ##
  ##                                                         ##
  ## ------------------------------------------------------- ##

  if(!is.null(fixed.cov)){
    Vd_fixed <- list()
    Vu_fixed <- list()
    for (i in 1:n.fixed.cov)
    {
      Vu_fixed[[i]] <- diag(me.fixed.cov[,i])
      Vd_fixed[[i]] <- diag(rep(var(fixed.cov[,i]), n)) - Vu_fixed[[i]]
      if(any(Vd_fixed[[i]] < 0)){
        Vd_fixed[[i]][Vd_fixed[[i]] < 0] <- 0
        warning("Vd contains negative variances, scaled up to 0. Does one or more predictor variables have a larger within-species measurement error than their among-species trait variance?")
      }
    }
  }else{
    Vd_fixed <- NULL
    Vu_fixed <- NULL
  }
  
  if(!is.null(random.cov)){
    Vu_random <- list()
    Vd_random <- lapply(sigma_squared, function(e) ta*e)
    for(i in seq_along(sigma_squared)){
      Vu_random[[i]] <- diag(me.random.cov[,i])
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
  
  list(sigma_squared = sigma_squared,
       theta.X = theta.X,
       Vu_given_x = Vu_given_x,
       Vu = Vu,
       Vx = Vx,
       Vd = Vd)
}