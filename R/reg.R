## Generalized inverse, or pseudoinverse
pseudoinverse <-
  function (m){
    msvd <- svd(m)
    if (length(msvd$d) == 0) {
      return(array(0, dim(m)[2:1]))
    }
    else {
      return(msvd$v %*% (1 / msvd$d * t(msvd$u)))
    }
  }

## Design matrix
slouch.modelmatrix <- function(a, hl, tree, observations, control, evolutionary= F){

  n <- length(tree$phy$tip.label)
  
  if(!is.null(observations$direct.cov) | !is.null(observations$random.cov)){
    if(evolutionary){
      rho <- 1
    }else{
      if(control$model == "ou"){
        rho <- (1 - (1 - exp(-a * tree$T.term))/(a * tree$T.term))
      }else{
        rho <- tree$T.term / 2
      }
    }

  }else{
    rho <- NULL
  }

  if(control$estimate.bXa){
    if (!is.null(observations$random.cov) | !is.null(observations$direct.cov)){
      bXa <- c(bXa = 1 - exp(-a * tree$T.term) - (1 - (1 - exp(-a * tree$T.term))/(a * tree$T.term)))
    }else{
      stop("bXa can not be estimated without continuous covariates.")
    }
  }else{
    bXa <- NULL
  }
  
  if(!is.null(observations$fixed.fact)){
    if(control$model == "ou"){
      c_regimes <- weight.matrix(tree$phy, a, tree$lineages)
    }else{
      c_regimes <- weight.matrix.brown(tree$lineages)
    }
    if(control$model == "bm" & control$estimate.Ya){
      X <- cbind(Ya = 1,
                 c_regimes)
    }else{
      X <- cbind(c_regimes, 
                 bXa)
    }
  }else{
    if(control$estimate.Ya){
      X <- cbind(Ya = exp(-a * tree$T.term),
                 b0 = 1-exp(-a * tree$T.term),
                 bXa)
    }else{
      X <- cbind("(Intercept)" = rep(1, length(tree$phy$tip.label)), 
                 bXa)
    }
    
  }
  
  if (!is.null(observations$fixed.fact) & control$interactions){
    ## Interactions
    direct.cov.list <- apply(observations$direct.cov, 2, list)
    interactions.list <- lapply(direct.cov.list, function(x) apply(c_regimes, 2, function(e) e * x[[1]]))
    interactions <- do.call(cbind, interactions.list)
    interactions.names.list <- lapply(colnames(observations$direct.cov), function(x) sapply(colnames(c_regimes), function(e) paste0(e, ":", x)))
    colnames(interactions) <- do.call(c, interactions.names.list)
    
    X <- cbind(X, 
               interactions,
               observations$random.cov*rho)
  }else{
    X <- cbind(X, 
               observations$direct.cov,
               observations$random.cov*rho)  
  }
 
  return(X)
}

#varcov_model <- function(hl, vy, a, beta1, which.direct.cov, which.random.cov, random.cov, T.term, direct.cov, Vu_given_x, mcov.random.cov, mcov.direct.cov, n, sigma_squared, phy, ta, tij, tja, mv.response){
varcov_model <- function(hl, sigma2_y, a, beta1, w_beta, tree, observations, seed, control){

  tij <- tree$tij
  tja <- tree$tja
  ta <- tree$ta
  n <- length(tree$phy$tip.label)
  
  if (hl == 0){
    
    Vt <- diag(rep(sigma2_y, times = n))
  }else if (a < 1e-14){
    if (!is.null(observations$random.cov)){
      s1 <- sum(seed$brownian_sigma_squared %*% ((beta1[w_beta$which.random.cov, ])^2))
      Vt <- sigma2_y * ta + s1 * ta * ((ta^2)/12 + tja*t(tja)/4)
    }else{
      Vt <- sigma2_y * ta
    }

  }else{
    if (!is.null(observations$random.cov)){
      s1 <- sum(seed$brownian_sigma_squared %*% ((beta1[w_beta$which.random.cov, ])^2))
      ti <- matrix(rep(tree$T.term, n), nrow=length(tree$phy$tip.label))
      
      term0 <- ((s1 + sigma2_y) / (2 * a)) * (1 - exp( -2 * a * ta)) * exp(-a * tij)
      term1 <- (1 - exp(-a * ti))/(a * ti) 
      term2 <- exp(-a * tja) * (1 - exp(-a * ti)) / (a * ti) 
      
      Vt <- term0 + s1*(ta * term1 * t(term1) - ((1 - exp(-a * ta))/a) * (term2 + t(term2)))
      
    }else{
      #Vt <- vy * (1 - exp(-2 * a * ta)) * exp(-a * tij)
      Vt <- sigma2_y/(2*a) * observations$closures$V_fixed_partial(a)
    }
  }
  return(Vt)
}

varcov_measurement <- function(observations, seed, beta1, hl, a, c_regimes, T.term, w_beta, control, evolutionary){

  if (evolutionary){
    rho2 <- 1
  }else{
    if(control$model == "ou"){
      rho2 <- (1 - (1 - exp(-a * T.term))/(a * T.term))^2
    }else{
      rho2 <- (T.term/2)^2
    }
  }

  
  
  if(length(w_beta$which.direct.cov) > 0 | length(w_beta$which.random.cov) > 0){
    ## Update measurement error in X.
    if (!is.null(observations$fixed.fact) & control$interactions){
      c_regime <- split(c_regimes, colnames(c_regimes)) #apply(c_regimes, 2, f)
      Vu_given_xi_direct <- head(seed$Vu_given_x, length(w_beta$which.direct.cov)/length(w_beta$which.regimes))
      Vu_given_xi_random <- tail(seed$Vu_given_x, length(w_beta$which.random.cov))
      Vu_given_xi <- c(do.call(c, lapply(Vu_given_xi_direct, function(e) lapply(w_beta$which.regimes, function(x) e))), 
                       Vu_given_xi_random)
      
    }else{
      c_regime <- 1
      Vu_given_xi <- seed$Vu_given_x
    }
    
    beta2_Vu_given_x_list <-  mapply(function(Vu_given_xi, b, w, rho_squared) {Vu_given_xi * b^2 * w * rho_squared} , 
                                     Vu_given_xi = Vu_given_xi, 
                                     b = beta1[c(w_beta$which.direct.cov, w_beta$which.random.cov),], 
                                     w = c(do.call(c, lapply(seq_along(colnames(observations$direct.cov)), function(e) c_regime)),
                                           lapply(seq_along(colnames(observations$random.cov)), function(e) 1)), ## interactions not supported for BM covariates
                                     rho_squared = c(lapply(seq_along(beta1[w_beta$which.direct.cov, ]), function(e) 1),
                                                     lapply(seq_along(beta1[w_beta$which.random.cov, ]), function(e) rho2)),
                                     SIMPLIFY = FALSE)
    beta2_Vu_given_x <- Reduce("+", beta2_Vu_given_x_list)
    
    # ## BROKEN!? Needs test.
    # calculate covariances between response and the stochastic predictor, to be subtracted in the diagonal of V
    if(sum(observations$mcov.random.cov) == 0){
      mcov.random <- 0
    }else{
      mcov.random <- rowSums(matrix(data=as.numeric(observations$mcov.random.cov)*t(kronecker(2*beta1[w_beta$which.random.cov,],
                                                                                        (1-(1-exp(-a*T.term))/(a*T.term)))), 
                             ncol = ncol(observations$random.cov)))
    }
    
    if(sum(observations$mcov.direct.cov) == 0){
      mcov.direct <- 0
    }else{
      mcov.direct <- rowSums(matrix(data=as.numeric(observations$mcov.direct.cov)*t(kronecker(2*beta1[w_beta$which.direct.cov,],
                                                                                             (1-(1-exp(-a*T.term))/(a*T.term)))), 
                                   ncol =  ncol(observations$direct.cov)))
    }
  }else{
    beta2_Vu_given_x <- 0
    mcov.random <- 0
    mcov.direct <- 0
  }
  
  V_me <- observations$mv.response + beta2_Vu_given_x - mcov.random - mcov.direct
  return(V_me)
}

reg <- function(par, tree, observations, control, seed, parameter_search = TRUE, gridsearch = FALSE){
  Y <- observations$Y

  if(is.numeric(par)){
    if(length(par) >= 2){
      par <- as.list(par)
    }else{
      par <- list(sigma2_y = par)
    }
  }
  
  if (control$model == "ou"){
    if (!is.null(par$hl)){
      hl <- par$hl  
      a <- log(2) / hl
    }else{
      a <- par$a
      hl <- log(2) / a
    }
    
    if(!is.null(par$vy)){
      vy <- par$vy
      sigma2_y <- vy*2*a 
    }else{
      sigma2_y <- par$sigma2_y
      vy <- sigma2_y/(2*a)
    }
    

  }else{
    #print(par)
    #sigma2_y <- if(length(par) == 1) par[[1]] else par$sigma2_y
    sigma2_y <- par$sigma2_y
    #sigma2_y <- unlist(par)
    hl <- Inf
    a <- 0
    vy <- NULL
  }
  
  
  n <- length(tree$phy$tip.label)
  
  if (control$verbose){
    cat(as.numeric(round(cbind(hl, (if(!is.null(vy)) vy else sigma2_y)), 4)))
    cat(" ")
  }
  
  X <- slouch.modelmatrix(a, hl, tree, observations, control, evolutionary = F)


  ## Initial OLS estimate of beta, disregarding variance covariance matrix
  beta1 <- matrix(stats::lm.fit(X, Y)$coefficients, ncol=1)
  
  n_direct <- length(observations$names.direct.cov)
  n_random <- length(observations$names.random.cov)
  regimes_all <- concat.factor(observations$fixed.fact, tree$regimes)
  n_regimes <- length(levels(regimes_all))
  
  
  #interaction = TRUE
  if (control$interactions){
    n_interaction <- n_direct * (n_regimes * control$interactions - 1)
  }else{
    n_interaction <- 0
  }
  
  beta1.descriptor <- c(rep("Intercept", length(beta1) - n_regimes - n_direct - n_random - n_interaction),
                        rep("regime", n_regimes),
                        rep("Instantaneous cov", n_direct + n_interaction),
                        rep("Random cov", n_random))
  
  which.intercepts <- which(beta1.descriptor == "Intercept")
  which.regimes <- which(beta1.descriptor == "regime")
  which.direct.cov <- which(beta1.descriptor == "Instantaneous cov")
  which.random.cov <- which(beta1.descriptor == "Random cov")

  w_beta <- list(which.intercepts = which.intercepts,
                 which.regimes = which.regimes,
                 which.random.cov = which.random.cov,
                 which.direct.cov = which.direct.cov)
  c_regimes <- X[ ,which.regimes]
  
  con.count <- 0
  repeat{
    Vt <- varcov_model(hl, sigma2_y, a, beta1, w_beta, tree, observations, seed, control)
    V_me <- varcov_measurement(observations, seed, beta1, hl, a, c_regimes, tree$T.term, w_beta, control, evolutionary = FALSE)
    V <- Vt + diag(V_me)
    
    L <- t(chol(V))

    # Note: Calculation of regression coefficients deviates from paper notation, but should be algebraically 
    # equivalent. This is done for an increase in computational speed, convenience of diagnosing problems
    # with better error-messages, and it is possibly more numerically stable.
    #
    # Instead of writing out the naive GLS beta estimator as-is, we have
    #
    # Linear transformation of both Y and model matrix, by the inverse of cholesky decomposition of V
    # With this method, the program will not crash even if X is singular. In this case, R's lm.fit() does not throw
    # error, but return coefficients that contain "NA", and subsequently the log-likelihood is evaluated as "NA".
  
    fit <- lm.fit(backsolve(L, X, upper.tri = FALSE), 
                  matrix(backsolve(L, Y, upper.tri = FALSE), ncol=1))
    beta.i <- matrix(fit$coefficients, ncol=1)
    
    ## Defensive & debug conditions
    if(any(is.na(beta.i[c(which.direct.cov, which.random.cov), ]))) {
      cat("\n")
      print(beta.i)
      warning("For hl = ", hl," and vy = ", vy," the gls estimate of beta contains \"NA\". Consider using different hl or vy.")
      print(as.numeric(round(cbind(hl, vy, NA, t(beta1)), 4)))
        return(list(support = NA,
                    par = par))
    }
      
    ## Check for convergence
    con.count <- con.count + 1
    if(all(abs(beta1[!is.na(beta1)] - beta.i[!is.na(beta.i)]) < control$convergence)) {
      beta1<-beta.i
      break
    }else{
      beta1<-beta.i
      
      if(con.count >= 50)
      {
        message("Warning, estimates did not converge after 50 iterations, last estimates printed out")
        break
      }
    }
  }

  
  log.det.V = 2*sum(log(diag(L))) ## Special case for positive definite matrix V, more numerically stable than (log*det(V)) routine for large V
  
  ## Log-likelihood
  sup1 <- - n/2*log(2*pi) - 0.5*log.det.V - 0.5*crossprod(fit$residuals)
  if(control$verbose){
    cat(as.numeric(round(cbind(sup1, t(beta1)), 4)))
    cat("\n")
  }

  if(parameter_search){
    return(list(support = sup1,
                par = par,
                gridsearch = gridsearch
                #V = V # This will cause memory overflow when n or Grid is large, or when less memory is available. O(n^2) + O(grid.vy * grid.hl)
                ))
  }else{
    
    ## Following is less optimized, more identical to publication syntax because it only has to run once.
    V.inverse <- solve(V)
    #beta1.var <- solve(crossprod(forwardsolve(L, X)))
    beta1.var <- pseudoinverse(crossprod(forwardsolve(L, X)))
    dimnames(beta1.var) <- list(colnames(X), colnames(X))
    
    if(control$model == "bm" & !is.null(observations$fixed.fact)){
      which.trend <- seq_along(levels(tree$regimes))
      trend_names <- colnames(X)[which.trend]
      trend_diff <- trend_diff_foo(tree, beta1, beta1.var, trend_names, which.trend)
    }else{
      trend_diff <- NULL
    }

    
    ## Evolutionary regression
    X0 <- slouch.modelmatrix(a = a, hl = hl, tree, observations, control, evolutionary = T)
    if(!is.null(observations$random.cov)){
      ev.beta1.var <- pseudoinverse(t(X0)%*%V.inverse%*%X0)
      dimnames(ev.beta1.var) <- list(colnames(X0), colnames(X0))
      
      ev.beta1 <- ev.beta1.var%*%(t(X0)%*%V.inverse%*%Y)
      
      
      beta_evolutionary <- c(list(coefficients = matrix(cbind(ev.beta1, sqrt(diag(ev.beta1.var))), 
                                             nrow=ncol(X0), 
                                             dimnames = list(colnames(X0), c("Predictions", "Std. error"))),
                       X = X0,
                       vcov = ev.beta1.var,
                       residuals = Y - (X0 %*% ev.beta1)),
                  bias_correction(ev.beta1, ev.beta1.var, Y, X0, V, w_beta, c_regimes, seed, observations, control))
    }else{
      ev.beta1.var <- NULL
      ev.beta1 <- NULL
      beta_evolutionary <- NULL
    }
    
    beta_primary <- c(list(coefficients = matrix(cbind(beta1, sqrt(diag(beta1.var))), 
                                            nrow=ncol(X), 
                                            dimnames = list(colnames(X), c("Estimates", "Std. error"))),
                      X = X,
                      residuals = Y - (X0 %*% beta1),
                      trend_diff = trend_diff,
                      vcov = beta1.var),
                 if (!is.null(observations$direct.cov) | !is.null(observations$random.cov)) bias_correction(beta1, beta1.var, Y, X, V, w_beta, c_regimes, seed, observations, control) else NULL)
    
    pred.mean <- X%*%beta1
    g.mean <- (t(rep(1, times = n)) %*% solve(V) %*% Y) / sum(solve(V))
    sst <- t(Y - c(g.mean)) %*% solve(V) %*% (Y - c(g.mean))
    sse <- t(Y - pred.mean) %*% solve(V) %*% (Y - pred.mean)
    r.squared <- (sst - sse) / sst

    return(list(support = sup1,
                par = par,
                V = V, 
                beta_primary = beta_primary,
                beta_evolutionary = beta_evolutionary,
                sst = sst,
                sse = sse,
                r.squared = r.squared,
                w_beta = w_beta))
  }
}


trend_diff_foo <- function(tree, beta1, beta1.var, trend_names, which.trend){
  if(length(which.trend) == 1){
    res <- cbind("Contrast" = "NA",
                 "Std. error" = "NA")
    return(res)
  }
  beta1 <- beta1[which.trend]
  beta1.var <- beta1.var[which.trend, which.trend]
  
  if(!is.null(tree$lineages[[1]]$lineage_regimes)){
    relevant <- upper.tri(beta1.var) | lower.tri(beta1.var)
    
    ## calculate which transitions are actually present on the tree
    transitions <- c()

    for (lineage in tree$lineages){
      anc <- tail(lineage$lineage_regimes, n = -1)
      desc <- head(lineage$lineage_regimes, n = -1)
      
      for (i in seq_along(anc)){
        if (anc[i] != desc[i]){
          trans <- paste(anc[i], "-", desc[i])
          transitions <- c(transitions, trans)
        }
      }
    }
    transitions <- unique(transitions)
  }else{
    transitions <- NULL
    relevant <- upper.tri(beta1.var)
  }
  
  
  se2 <- diag(beta1.var)
  contrast <- sapply(beta1, function(e) e - beta1)
  v <- sapply(se2, function(e) e + se2)
  names_matrix <- sapply(trend_names, function(e) paste(e, "-", trend_names))
  
  se_contrast <- sqrt(v[relevant] -2*beta1.var[relevant])
  res <- cbind("Contrast" = contrast[relevant],
               "Std. error" = se_contrast)
  rownames(res) <- names_matrix[relevant]
  
  if(!is.null(transitions)){
    ## Filter for only the transitions that are on the tree
    res <- res[rownames(res) %in% transitions, ]  
  }

  
  
  return(res)
}

is.diag <- function(m){
  all(m[lower.tri(m)] == 0, m[upper.tri(m)] == 0)
}