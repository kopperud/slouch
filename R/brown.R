## Extra functions

sigma.X.estimate <-
  function (phy, ta, predictor, mv.predictor) {
    predictor <- matrix(predictor, nrow = length(phy$tip.label))
    
    N <- length(phy$tip.label)
    v <- ta # Time from root to most recent ancestor
    w <- matrix(data = 1, nrow = N, ncol = 1)
    me <- diag(mv.predictor)
    dat <- predictor
    beta <- solve(t(w) %*% solve(v) %*% w) %*% (t(w) %*% solve(v) %*% dat)
    e <- dat - c(beta)
    sigma_squared <- as.numeric((t(e) %*% solve(v) %*% e) / (N-1))
    repeat{
      beta <- solve(t(w) %*% solve(v + me/sigma_squared) %*% w) %*% (t(w) %*% solve(v + me/sigma_squared) %*% dat)
      e <- dat - c(beta)
      theta <- beta
      sigma_squared1 <- (t(e) %*% solve(v + me/sigma_squared) %*% e) / (N-1)
      if (abs(as.numeric(sigma_squared1) - sigma_squared) <= 0.0000001 * sigma_squared){
        break
      } 
      sigma_squared <- as.numeric(sigma_squared1)
    }
    return(list(mean = as.numeric(theta),
                sigma_squared = as.numeric(sigma_squared)))
  }