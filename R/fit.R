## Internal model fitting function
.slouch.fit <- function(phy,
                        species,
                        hl_values, 
                        a_values,
                        vy_values, 
                        sigma2_y_values,
                        response, 
                        mv.response, 
                        fixed.fact,
                        direct.cov, 
                        mv.direct.cov, 
                        mcov.direct.cov, 
                        random.cov, 
                        mv.random.cov, 
                        mcov.random.cov,
                        estimate.Ya,
                        estimate.bXa,
                        interactions,
                        hessian,
                        model,
                        support, 
                        convergence,
                        nCores,
                        hillclimb,
                        lower,
                        upper,
                        verbose,
                        names.direct.cov,
                        names.random.cov,
                        name.response)
{
  if(is.null(species)){
    stop("Use argument \"species\" to make sure the order of the data correctly lines up with the tree. See example.")
  }
  if(!all(species == phy$tip.label)){
    stop("Argument \"species\" must have the exact same species names as phy$tip.label, and in the same order.")
  }
  
  ## Checks, defensive conditions
  stopifnot(ape::is.rooted(phy))
  
  if(estimate.Ya | estimate.bXa){
    if(ape::is.ultrametric(phy)){
      warning("Your tree looks ultrametric - ancestral state coefficients may cause the model matrix to be singular. Consider not using \"estimate.Ya\" and\\or \"estimate.bXa\".")
    }
  }
  
  # SET DEFAULTS IF NOT SPECIFIED
  if(is.null(mv.response)){
    mv.response <- rep(0, times= length(response))
  }
  
  # SPECIFY COMPONENTS THAT ARE COMMON TO ALL MODELS
  
  Y <- response
  n <- length(Y)
  
  if(!is.null(fixed.fact)){
    if(is.null(phy$node.label)){
      stop("For categorical variables, the regimes corresponding to their primary optima need to be painted on all of the branches in the tree, and assigned to phy$node.label - use plot(phy) & nodelabels(phy$node.label) to see whether they are correct. See example.")
    }
    
    regimes_internal <- phy$node.label
    if(estimate.Ya & model == "ou"){
      tmp <- as.character(regimes_internal)
      tmp[1] <- "Ya"
      regimes_internal <- factor(tmp)
    }
    regimes_tip <- fixed.fact
    
    regimes <- concat.factor(regimes_tip, regimes_internal)
    lineages <- lapply(1:n, function(e) lineage.constructor(phy, e, regimes)) #; names(lineages) <- phy$tip.label
  }else{
    regimes <- lineages <- NULL
  }
  
  mrca1 <- ape::mrca(phy)
  times <- ape::node.depth.edgelength(phy)
  ta <- matrix(times[mrca1], nrow=n, dimnames = list(phy$tip.label, phy$tip.label))
  T.term <- times[1:n]
  tia <- times[1:n] - ta
  tja <- t(tia)
  tij <- tja + tia
  
  ############################################################################
  
  ##          Make sure variable matrices are of correct dimensions         ##
  
  ############################################################################
  
  if(!is.null(direct.cov)){
    if(ncol(as.matrix(direct.cov))==1) {
      #names.direct.cov <- deparse(substitute(direct.cov))
      direct.cov <- matrix(direct.cov, nrow = length(phy$tip.label), dimnames = list(NULL, names.direct.cov))
      
    }else{
      direct.cov <- as.matrix(direct.cov)
      stopifnot(!is.null(colnames(direct.cov)))
      #names.direct.cov <- colnames(direct.cov)
    }
    if(is.null(mv.direct.cov)){
      mv.direct.cov <- matrix(0, nrow = n, ncol = ncol(direct.cov))
    }else{
      mv.direct.cov <- cbind(mv.direct.cov)
    }
    
    if(is.null(mcov.direct.cov)){
      mcov.direct.cov <- matrix(data = 0, nrow = n, ncol = ncol(direct.cov))
    }else{
      mcov.direct.cov <- as.matrix(mcov.direct.cov)
    }
    
  }else{
    names.direct.cov <- NULL
  }
  
  if(!is.null(random.cov)){
    if(ncol(as.matrix(random.cov)) == 1) {
      #names.random.cov <- deparse(substitute(random.cov))
      random.cov <- matrix(random.cov, nrow = length(phy$tip.label), dimnames = list(NULL, names.random.cov))
    }
    else{
      random.cov <- as.matrix(random.cov)
      stopifnot(!is.null(colnames(random.cov)))
      #names.random.cov <- colnames(random.cov)
    }
    if(is.null(mv.random.cov)){
      mv.random.cov <- matrix(0, nrow = n, ncol = ncol(random.cov))
    }else{
      mv.random.cov <- cbind(mv.random.cov)
    }
    
    if(is.null(mcov.random.cov)){
      mcov.random.cov <- matrix(data = 0, nrow = n, ncol = ncol(random.cov))
    }else{
      mcov.random.cov <- as.matrix(mcov.random.cov)
    }
  }else{
    names.random.cov <- NULL
  }
  
  
  ## Cluster all parameters concerning phylogenetic tree
  tree <- list(phy = phy,
               T.term = T.term,
               tia = tia,
               tja = tja,
               ta = ta,
               tij = tij,
               times = times,
               lineages = lineages,
               regimes = regimes)
  
  ## Cluster parameters concerning the type of model being run
  observations <- list(response = response,
                       mv.response = mv.response,
                       fixed.fact = fixed.fact,
                       direct.cov = direct.cov,
                       mcov.direct.cov = mcov.direct.cov,
                       random.cov = random.cov,
                       mcov.random.cov = mcov.random.cov,
                       Y = Y,
                       names.direct.cov = names.direct.cov,
                       names.random.cov = names.random.cov,
                       closures = list(V_fixed_partial = memoise::memoise(function(a) (1 - exp(-2 * a * ta)) * exp(-a * tij))))
  
  control <- list(verbose = verbose,
                  estimate.Ya = estimate.Ya,
                  estimate.bXa = estimate.bXa,
                  support = support,
                  convergence = convergence,
                  model = model,
                  name.response = name.response,
                  interactions = interactions)
  
  
  seed <- seed(phy, ta, direct.cov, mv.direct.cov, random.cov, mv.random.cov)
  coef.names <- colnames(slouch.modelmatrix(a = 1, hl = 1, tree, observations, control, evolutionary = T))
  
  ## Uniform random start values for hl and vy, in case hillclimber is used
  if(model == "bm"){
    if(is.null(sigma2_y_values)){
      sigma2_y_values <- stats::runif(1, 0, stats::var(response))
    }
  }else{
    ## Random numbers
    if (is.null(vy_values) & is.null(sigma2_y_values)){
      X0 <- slouch.modelmatrix(1000000, 0, tree, observations, control, evolutionary= F)
      model0  <- stats::lm.fit(X0, observations$response)
      vy_values <- stats::var(model0$residuals)
      #vy_values <- stats::runif(1, 0, stats::var(response))
    }else{
      
    }
    if (is.null(hl_values) & is.null(a_values)){
      hl_values <- stats::runif(1, 0, max(times))
    }
  }
  
  
  if(verbose){
    message("GRID SEARCH PARAMETER SUPPORT")
    cat(c("     hl     ", if(model == "ou") "vy    " else "sigma2_y    ", "support", c(coef.names), "\n"))
  }
  
  #############
  if(model == "bm"){
    gridpar <- lapply(sigma2_y_values, function(e) list(sigma2_y = e
    ))
  }else{
    gridlist <- list(hl = hl_values,
                     a = a_values,
                     vy = vy_values,
                     sigma2_y = sigma2_y_values)
    gridlist[sapply(gridlist, is.null)] <- NULL ## Remove null entries
    gridm <- expand.grid(gridlist)
    gridpar <- apply(gridm, 1, function(e) as.list(e))
  }
  
  time0 <- Sys.time()
  
  if(nCores > 1 & nCores %% 1 == 0){
    if(!"package:parallel" %in% search()){
      stop("For multicore, please load package with library(parallel)")
    }
    if(.Platform$OS.type == "unix"){
      grid <- parallel::mclapply(gridpar,
                                 function(e) reg(e, tree, observations, control, seed, gridsearch = TRUE),
                                 mc.cleanup = TRUE,
                                 mc.cores = nCores)
    }else{
      cl <- parallel::makeCluster(getOption("cl.cores", nCores))
      parallel::setDefaultCluster(cl)
      parallel::clusterExport(cl, c("tree", "observations", "control", "seed"), envir = environment())
      grid <- parallel::parLapply(cl, gridpar, function(e) reg(e, tree, observations, control, seed, gridsearch = TRUE))
      parallel::stopCluster(cl)
    }
  }else{
    grid <- lapply(gridpar, reg, tree, observations, control, seed, gridsearch = TRUE)
  }
  
  sup2_grid <- sapply(grid, function(e) e$support)
  which1 <- which.max(sup2_grid)
  ml_grid <- grid[[which1]]$support
  
  if(hillclimb){
    
    if(verbose){
      Sys.sleep(0.2)
      cat("\n")
      print(paste("Start hillclimb parameter estimation routine, method", if(model == "ou") "L-BFGS-B." else "Brent."))
      cat("\n")
    }
    
    hcenv <- environment()
    hcenv$k <- 0
    climblog <- list()
    par <- c(grid[[which1]]$par)
    
    if(model == "ou"){
      parscale <- c(max(T.term), var(response))
    }else{
      parscale <- var(response)
    }
    
    optimout <- stats::optim(
      par = par,
      fn = function(e, ...){hcenv$k <- hcenv$k +1; tmp <- reg(e, tree, observations, control, seed, ...); hcenv$climblog[[toString(hcenv$k)]] <- tmp; return(tmp$support) }, ## Ugly environment hack to log the hillclimber. Impure function
      parameter_search = TRUE,
      lower = lower,
      upper = upper,
      method = if(model == "ou") "L-BFGS-B" else "Brent",
      control = list(parscale = parscale,
                     fnscale = -1),
      hessian = hessian
    )

    
    climblog2 <- sapply(names(climblog[[1]]$par), 
                        function(x) sapply(climblog, 
                                           function(e) e$par[[x]]),
                        USE.NAMES = TRUE,
                        simplify = FALSE)
    climblog2[sapply(climblog2, is.null)] <- NULL
    
    ## Matrix for plotting the route of hillclimber
    climblog_df <- data.frame(index = seq_along(climblog),
                              loglik = sapply(climblog, function(e) e$support),
                              climblog2)
    
  }else{
    optimout <- NULL
    climblog_df <- NULL
    climblog <- NULL
  }
  parameter_space <- c(grid, climblog)
  sup2 <- sapply(parameter_space, function(e) e$support)
  ml <- max(sup2, na.rm = T)
  
  
  ## Find the regression for which the support value is maximized
  par = parameter_space[[which.max(sup2)]]$par
  
  ## Repeat regression at a, vy for which logLik is maximized
  fit <- reg(par, tree, observations, control, seed, parameter_search = FALSE)
  
  if(verbose){
    print(paste0("Parameter search done after ",round((Sys.time() - time0), 3)," s."))
  }
  
  ############################
  
  if (model == "ou"){
    alpha = log(2) / fit$par$hl
    evolpar <- par
  }else{
    evolpar <- list(
      sigma2_y <- fit$par$sigma2_y
    )
    names(evolpar) <- "sigma2_y"
  }
  
  
  if(!is.null(random.cov)){
    brownian_predictors <- matrix(data = rbind(seed$brownian_mean, seed$brownian_sigma_squared), 
                                  byrow = TRUE,
                                  ncol = 2, 
                                  dimnames = list(names.random.cov,
                                                  c("Phylogenetic mean", "Diffusion variance"))
    )
  }else{
    brownian_predictors <- NULL
  }
  
  n.par <- length(coef.names) + if (model == "ou") 2 else 1
  
  aic = -2*ml+2*n.par
  aicc = aic + (2*n.par*(n.par+1))/(n-n.par-1)
  sic = -2*ml+log(n)*n.par
  r.squared = fit$r.squared*100
  
  modfit <- list(Support = ml,
                 AIC = -2*ml+2*n.par,
                 AICc = aic + (2*n.par*(n.par+1))/(n-n.par-1),
                 SIC = -2*ml+log(n)*n.par,
                 "R squared" = fit$r.squared,
                 SST = fit$sst,
                 SSE = fit$sse,
                 "N (params)" = n.par)
  
  if(model == "ou"){
    if(!is.null(c(hl_values, a_values)) & !is.null(c(vy_values, sigma2_y_values))){
      if(length(c(hl_values, a_values)) > 1 | length(c(vy_values, sigma2_y_values)) > 1){
        
        z <- matrix(sapply(grid, function(e) e$support), nrow=length(c(hl_values, a_values)), byrow=F)
        if(!(all(is.na(z) | is.infinite(z)))){
          z <- z - ml
          z[abs(z) >= support] <- -2
          
          supportplot <- list(a = a_values,
                              hl = hl_values,
                              vy = vy_values,
                              sigma2_y = sigma2_y_values,
                              z = z)
        }else{
          warning("All support values in grid either NA or +/-Inf - Can't plot.")
          supportplot <- NULL
        }
      }else{
        supportplot <- NULL
      }
    }else{
      supportplot <- NULL
    }
    
  }else{
    if(length(sigma2_y_values) > 1){
      supportplot <- data.frame(sigma2_y = sapply(grid, function(e) e$par$sigma2_y),
                                loglik = sapply(grid, function(e) e$support))
      
    }else{
      supportplot <- NULL
    }
  }
  supported_range <- support_interval(grid, ml, support)
  
  
  result <- list(parameter_space = parameter_space,
                 tree = tree,
                 modfit = modfit,
                 supportplot = supportplot,
                 climblog_df = climblog_df,
                 brownian_predictors = brownian_predictors,
                 beta_primary = fit$beta_primary,
                 beta_evolutionary = fit$beta_evolutionary,
                 evolpar = evolpar,
                 supported_range = supported_range,
                 n.par = n.par,
                 V = fit$V,
                 fixed.fact = fixed.fact,
                 control = control,
                 hessian = optimout$hessian,
                 w_beta = fit$w_beta)
  class(result) <- c("slouch", class(result))
  return(result)
  
}

support_interval <- function(grid, ml, support){
  loglik <- sapply(grid, function(x) x$support)
  supported <- grid[loglik >= ml - support]
  
  if (length(supported) > 0){
    supported2 <- sapply(names(grid[[1]]$par), 
                         function(x) sapply(supported, 
                                            function(e) e$par[[x]]), 
                         USE.NAMES = TRUE, 
                         simplify = F)
    supported2[is.null(supported2)] <- NULL
    supported_range <- t(sapply(supported2, range))
    colnames(supported_range) <- c("Minimum", "Maximum")
  }else{
    supported_range <- NULL
  }
  return(supported_range)
}