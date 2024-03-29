## ------------------------------------------ ##
##              This is the worst             ##
##         part of the code. 27 aug 2016      ##
## ------------------------------------------ ##


seed <- function(phy, ta, direct.cov, mv.direct.cov, random.cov, mv.random.cov){
  n <- length(phy$tip.label)
  
  if(!is.null(random.cov)){
    brownian <- list()
    for (i in 1:ncol(random.cov)){
      y <- random.cov[,i]
      y_me <- mv.random.cov[,i] 
      
      brownian[[i]] <-  sigma.X.estimate(phy, ta, y, y_me)
    }
    names(brownian) <- colnames(random.cov)

    sigma_squared <- sapply(brownian, function(x) x$sigma_squared)
    brownian_mean <- sapply(brownian, function(x) x$mean)
  }else{
    sigma_squared <- NULL
    brownian_mean <- NULL
  }

  ## ------------------------------------------------------- ##
  ##                                                         ##
  ##                Defining Vu, Vd, Vx, Vu|x                ##
  ##     See equation (10) in Hansen & Bartoszek 2012        ##
  ##                                                         ##
  ## ------------------------------------------------------- ##

  if(!is.null(direct.cov)){
    Vd_fixed <- list()
    Vu_fixed <- list()
    for (i in 1:ncol(direct.cov))
    {
      Vu_fixed[[i]] <- diag(mv.direct.cov[,i])
      Vd_fixed[[i]] <- diag(rep(stats::var(direct.cov[,i]), n)) - Vu_fixed[[i]]
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
      Vu_random[[i]] <- diag(mv.random.cov[,i])
    }
  }else{
    Vd_random <- NULL
    Vu_random <- NULL
  }
  
  Vu <- c(Vu_fixed, Vu_random)
  Vd <- c(Vd_fixed, Vd_random)
  
  if(!is.null(random.cov) | !is.null(direct.cov)){
    Vu_given_x <- list()
    Vx <- list()
    for (i in 1:length(Vd)){
      Vx[[i]] <- Vd[[i]] + Vu[[i]]
      Vu_given_x[[i]] <- diag(Vu[[i]] - (Vu[[i]] %*% solve(Vx[[i]]) %*% Vu[[i]])) ## Turn diagonal matrix into vector
    }
  }else{
    Vu_given_x = NULL
    Vx <- NULL
  }
  
  list(brownian_sigma_squared = sigma_squared,
       brownian_mean = brownian_mean,
       Vu_given_x = Vu_given_x,
       Vu = Vu,
       Vx = Vx,
       Vd = Vd)
}