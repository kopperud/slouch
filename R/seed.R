## ------------------------------------------ ##
##              This is the worst             ##
##         part of the code. 27 aug 2016      ##
## ------------------------------------------ ##


seed <- function(phy, ta, fixed.cov, me.fixed.cov, random.cov, me.random.cov){
  n <- length(phy$tip.label)
  
  if(!is.null(random.cov)){
    brownian <- mapply(function(y, y_me) sigma.X.estimate(phy, ta, y, y_me),
                       y = split(t(random.cov), colnames(random.cov)),
                       y_me = split(t(me.random.cov), colnames(random.cov)),
                       SIMPLIFY = FALSE)
    
    sigma_squared <- sapply(brownian, function(x) x$sigma_squared)
    theta.X <- sapply(brownian, function(x) x$mean)
  }else{
    sigma_squared <- NULL
    theta.X <- NULL
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
    for (i in 1:ncol(fixed.cov))
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