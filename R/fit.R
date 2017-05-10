#' Function to fit Ornstein-Uhlenbeck models
#'
#' @param phy an object of class 'phylo', must be rooted.
#' @param species a character vector of species tip labels, typically the "species" column in a data frame. This column needs to be an exact match and same order as phy$tip.label
#' @param hl_values a vector of candidate phylogenetic half-life values to be evaluated in grid search. Optional.
#' @param vy_values a vector of candidate stationary variances for the response trait, to be evaluated in grid search. Optional.
#' @param response a numeric vector of a trait to be treated as response variable
#' @param me.response numeric vector of the observational variances of each response trait. E.g if response is a mean trait value, me.response is the within-species SE of the mean.
#' @param fixed.fact factor of regimes on the terminal edges of the tree, in same order as species. If this is used, phy$node.label needs to be filled with the corresponding internal node regimes, in the order of node indices (root: n+1),(n+2),(n+3), ...
#' @param fixed.cov Direct effect independent variables
#' @param me.fixed.cov Observational variances for direct effect independent variables. Must be the same shape as fixed.cov
#' @param mecov.fixed.cov .
#' @param random.cov Independent variables each modeled as a brownian motion
#' @param me.random.cov Observational variances for the brownian covariates. Must be the same shape as random.cov
#' @param mecov.random.cov .
#' @param estimate.Ya a logical value indicathing whether "Ya" should be estimated. If true, the intercept K = 1 is expanded to Ya = exp(-a*T.term) and b0 = 1-exp(-a*T.term). If models with categorical covariates are used, this will instead estimate a separate primary optimum for the root niche, "Ya". This only makes sense for non-ultrametric trees. If the tree is ultrametric, the model matrix becomes singular.
#' @param estimate.bXa a logical value indicathing whether "bXa" should be estimated. If true, bXa = 1-exp(-a*T.term) - (1-(1-exp(-a*T.term))/(a*T.term)) is added to the model matrix, estimating b*Xa. Same requirements as for estimating Ya.
#' @param support a scalar indicating the size of the support set, defaults to 2 units of log-likelihood.
#' @param convergence threshold of iterative GLS estimation for when beta is considered to be converged.
#' @param nCores number of CPU cores used in grid-search. If 2 or more cores are used, all print statements are silenced during grid search. If performance is critical it is recommended to compile and link R to a multithreaded BLAS, since most of the heavy computations are common matrix operations. Even if a singlethreaded BLAS is used, this may or may not improve performance, and performance may vary with OS.
#' @param hillclimb logical, whether to use hillclimb parameter estimation routine or not. This routine (L-BFGS-B from optim()) may be combined with the grid-search, in which case it will on default start on the sigma and halflife for the local ML found by the grid-search.
#' @param hillclimb_start numeric vector of length 2, c(hl, vy), to specify where the hillclimber routine starts.
#' @param lower lower bounds for the optimization routine, defaults to c(0,0). First entry in vector is half-life, second is stationary variance. When running direct effect models without observational error, it may be useful to specify a positive lower bounds for the stationary variance, e.g c(0, 0.001), since the residual variance-covariance matrix is degenerate when sigma = 0.
#' @param upper upper bounds for the optimization routine, defaults to positive infinite.
#' @param verbose a logical value indicating whether to print a summary in each iteration of parameter search. May be useful when diagnosing unexpected behaviour or crashes.
#'
#' @return An object of class 'slouch'
#' @export
slouch.fit<-function(phy,
                     species = NULL,
                     hl_values = NULL, 
                     vy_values = NULL, 
                     response, 
                     me.response=NULL, 
                     fixed.fact=NULL,
                     fixed.cov=NULL, 
                     me.fixed.cov=NULL, 
                     mecov.fixed.cov=NULL, 
                     random.cov=NULL, 
                     me.random.cov=NULL, 
                     mecov.random.cov=NULL,
                     estimate.Ya = FALSE,
                     estimate.bXa = FALSE,
                     support = 2, 
                     convergence = 0.000001,
                     nCores = 1,
                     hillclimb = FALSE,
                     hillclimb_start = NULL,
                     lower = c(0,0),
                     upper = c(Inf, Inf),
                     verbose = FALSE)
{
  if(is.null(species)){
    stop("Use argument \"species\" to make sure the order of the data correctly lines up with the tree. See example.")
  }
  if(!all(species == phy$tip.label)){
    stop("Argument \"species\" must have the exact same species names as phy$tip.label, and in the same order.")
  }
  
  ## Checks, defensive conditions
  #stopifnot(intercept == "root" | is.null(intercept))
  stopifnot(ape::is.rooted(phy))
  stopifnot((is.numeric(hillclimb_start) & length(hillclimb_start) == 2) | is.null(hillclimb_start))
  if((is.null(hl_values) | is.null(vy_values)) & !hillclimb){
    stop("Choose at minimum a 1x1 grid, or use the hillclimber routine.")
  }
  
  # SET DEFAULTS IF NOT SPECIFIED
  if(is.null(me.response)){
    me.response<-diag(rep(0, times= length(response)))
  }else{
    me.response<-diag(me.response)
  }
  
  # SPECIFY COMPONENTS THAT ARE COMMON TO ALL MODELS
  
  Y <- response
  n <- length(Y)
  
  if(!is.null(fixed.fact)){
    if(is.null(phy$node.label)){
      stop("For categorical variables, the regimes corresponding to their primary optima need to be painted on all of the branches in the tree, and assigned to phy$node.label - use plot(phy) & nodelabels(phy$node.label) to see whether they are correct. See example.")
    }
    
    regimes_internal <- phy$node.label
    if(estimate.Ya){
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
  
  h.lives<-matrix(data=0, nrow=length(hl_values), ncol=length(vy_values))
  hl_values<-rev(hl_values)
  
  ############################################################################
  
  ##          Make sure variable matrices are of correct dimensions         ##
  
  ############################################################################
  
  if(!is.null(fixed.cov)){
    if(ncol(as.matrix(fixed.cov))==1) {
      names.fixed.cov <- deparse(substitute(fixed.cov))
      fixed.cov <- matrix(fixed.cov, nrow = length(phy$tip.label), dimnames = list(NULL, names.fixed.cov))

    }else{
      fixed.cov <- as.matrix(fixed.cov)
      stopifnot(!is.null(colnames(fixed.cov)))
      names.fixed.cov <- colnames(fixed.cov)
    }
    if(is.null(me.fixed.cov)){
      me.fixed.cov <- matrix(0, nrow = n, ncol = ncol(fixed.cov))
    }else{
      me.fixed.cov <- cbind(me.fixed.cov)
    }
    
    if(is.null(mecov.fixed.cov)){
      mecov.fixed.cov <- matrix(data = 0, nrow = n, ncol = ncol(fixed.cov))
    }else{
      mecov.fixed.cov <- as.matrix(mecov.fixed.cov)
    }
    
  }else{
    names.fixed.cov <- NULL
  }
  
  if(!is.null(random.cov)){
    if(ncol(as.matrix(random.cov)) == 1) {
      names.random.cov <- deparse(substitute(random.cov))
      random.cov <- matrix(random.cov, nrow = length(phy$tip.label), dimnames = list(NULL, names.random.cov))
    }
    else{
      random.cov <- as.matrix(random.cov)
      stopifnot(!is.null(colnames(random.cov)))
      names.random.cov <- colnames(random.cov)
    }
    if(is.null(me.random.cov)){
      me.random.cov <- matrix(0, nrow = n, ncol = ncol(random.cov))
    }else{
      me.random.cov <- cbind(me.random.cov)
    }
    
    if(is.null(mecov.random.cov)){
      mecov.random.cov <- matrix(data = 0, nrow = n, ncol = ncol(random.cov))
    }else{
      mecov.random.cov <- as.matrix(mecov.random.cov)
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
  pars <- list(response = response,
               me.response = me.response,
               fixed.fact = fixed.fact,
               fixed.cov = fixed.cov,
               mecov.fixed.cov = mecov.fixed.cov,
               random.cov = random.cov,
               mecov.random.cov = mecov.random.cov,
               Y = Y,
               names.fixed.cov = names.fixed.cov,
               names.random.cov = names.random.cov)
  
  control <- list(verbose = verbose,
                  estimate.Ya = estimate.Ya,
                  estimate.bXa = estimate.bXa,
                  support = support,
                  convergence = convergence)
  
  
  seed <- seed(phy, ta, fixed.cov, me.fixed.cov, random.cov, me.random.cov)
  coef.names <- colnames(slouch.modelmatrix(a = 1, hl = 1, tree, pars, control, is.opt.reg = TRUE))
  
  ## Uniform random start values for hl and vy, in case hillclimber is used
  if (is.null(hl_values)){
    hl_values <- stats::runif(1, 0, max(times))
  }
  if (is.null(vy_values)){
    vy_values <- stats::runif(1, 0, stats::var(response))
  }
  
  if(verbose){
    message("GRID SEARCH PARAMETER SUPPORT")
    cat(c("     hl     ", "vy    ", "support", c(coef.names), "\n"))
  }
  
  #############
  vector_hl_vy <- cbind(sort(rep(hl_values, length(vy_values)), decreasing = TRUE), rep(vy_values, length(hl_values)))
  time0 <- Sys.time()
  
  if(nCores > 1 & nCores %% 1 == 0){
    if(!"package:parallel" %in% search()){
      stop("For multicore, please load package with library(parallel)")
    }
    if(.Platform$OS.type == "unix"){
      list_hl_vy <- lapply(seq_len(nrow(vector_hl_vy)), function(e) vector_hl_vy[e,])
      grid_support <- parallel::mclapply(list_hl_vy, 
                                         function(e) reg(e, tree, pars, control, seed),
                                         mc.cleanup = TRUE,
                                         mc.cores = nCores)
    }else{
      cl <- parallel::makeCluster(getOption("cl.cores", nCores))
      parallel::setDefaultCluster(cl)
      parallel::clusterExport(cl, c("tree", "pars", "control", "seed"), envir = environment())
      grid_support <- parallel::parApply(cl, vector_hl_vy, 1, function(e) reg(e, tree, pars, control, seed))
      parallel::stopCluster(cl)
    }
  }else{
    grid_support <- apply(vector_hl_vy, 1, reg, tree, pars, control, seed)
  }
  
  sup2_grid <- sapply(grid_support, function(e) e$support)
  ml_grid <- max(stats::na.exclude(sup2_grid))
  
  if(hillclimb){
    if(verbose){
      Sys.sleep(0.2)
      cat("\n")
      print("Start hillclimb parameter estimation routine, method L-BFGS-B")
      cat("\n")
    }
    
    hcenv <- environment()
    hcenv$k <- 0
    climblog <- list()
    if(is.null(hillclimb_start)){
      hillclimb_start <- vector_hl_vy[which.max(sup2_grid),]
    }
    hl_vy_est <- stats::optim(
      par = hillclimb_start,
      fn = function(e, ...){hcenv$k <- hcenv$k +1; tmp <- reg(e, tree, pars, control, seed, ...); hcenv$climblog[[toString(hcenv$k)]] <- tmp; return(tmp$support) }, ## Ugly environment hack to log the hillclimber. Impure function
      gridsearch = TRUE,
      lower = lower,
      upper = upper,
      method = "L-BFGS-B",
      control = list(parscale = c(max(T.term), var(response)),
                     fnscale = -0.1)
    )
    
    ## Matrix for plotting the route of hillclimber
    climblog_matrix <- data.frame(index = 1:length(climblog), 
                                  hl = sapply(climblog, function(e) e$hl_vy[[1]]), 
                                  vy = sapply(climblog, function(e) e$hl_vy[2]), 
                                  loglik = sapply(climblog, function(e) e$support))
  }else{
    climblog_matrix <- NULL
    climblog <- NULL
  }
  parameter_space <- c(grid_support, climblog)
  sup2 <- sapply(parameter_space, function(e) e$support)
  ml <- max(na.exclude(sup2))
  
  gof <- matrix(sup2_grid, ncol=length(vy_values), byrow=TRUE, dimnames = list(hl_values, vy_values))
  gof <- ifelse(gof <= ml-support, ml-support, gof) - ml
  
  
  ## Find the regression for which the support value is maximized
  besthl_vy = parameter_space[[which.max(sup2)]]$hl_vy
  
  ## Repeat regression at a, vy for which logLik is maximized
  fit <- reg(besthl_vy, tree, pars, control, seed, gridsearch=FALSE)
  
  if(verbose){
    print(paste0("Parameter search done after ",round((Sys.time() - time0), 3)," s."))
  }
  
  ############################
  
  alpha <- log(2) / fit$hl_vy[1]
  
  oupar <- matrix(c(alpha, 
                    fit$hl_vy[1], 
                    fit$hl_vy[2], 
                    mean((1-(1-exp(-alpha*T.term))/(alpha*T.term)))), 
                  ncol=1, 
                  dimnames=list(c("Rate of adaptation", "Phylogenetic half-life", "Stationary variance", "Phylogenetic correction factor"), "Estimate"))
  
  if(!is.null(random.cov)){
    brownian_predictors <- matrix(data = rbind(seed$theta.X, seed$sigma_squared), 
                                  nrow = 2, 
                                  ncol = ncol(random.cov), 
                                  dimnames = list(c("Predictor theta", "Predictor variance"), 
                                                  names.random.cov)
                                  )
  }else{
    brownian_predictors <- NULL
  }
  
  n.par <- length(coef.names) + 2
  modfit<-matrix(data=0, nrow=7, ncol=1, dimnames=list(c("Support", "AIC", "AICc", "SIC", "r squared", "SST", "SSE"),("Value")))
  modfit[1,1] = ml
  modfit[2,1] = -2*ml+2*n.par
  modfit[3,1] = modfit[2,1]+(2*n.par*(n.par+1))/(n-n.par-1)
  modfit[4,1] = -2*ml+log(n)*n.par
  modfit[5,1] = fit$r.squared*100
  modfit[6,1] = fit$sst
  modfit[7,1] = fit$sse
  
  
  if(length(hl_values) > 1 && length(vy_values) > 1){
    ## All hl + vy in the support interval
    hlsupport <- ifelse(sup2_grid <= ml - support, NA, sapply(grid_support, function(e) e$hl_vy[1]))
    vysupport <- ifelse(sup2_grid <= ml - support, NA, sapply(grid_support, function(e) e$hl_vy[2]))
    
    hlvy_grid_interval <- matrix(c(min(hlsupport, na.rm = TRUE), min(vysupport, na.rm = TRUE),
                                   max(hlsupport, na.rm = TRUE), max(vysupport, na.rm = TRUE)),
                                 ncol = 2, nrow = 2,
                                 dimnames = list(c("Phylogenetic half-life", "Stationary variance"), c("Minimum", "Maximum")))
    
    if(!(all(is.na(gof) | is.infinite(gof)))){
      h.lives <- matrix(0, nrow=length(hl_values), ncol=length(vy_values), dimnames = list(rev(hl_values), vy_values))
      for(i in 1:length(vy_values)){
        h.lives[,i]=rev(gof[,i])
      }
      
      supportplot = list(hl = rev(hl_values),
                         vy = vy_values,
                         z = h.lives)
    }else{
      warning("All support values in grid either NA or +/-Inf - Can't plot.")
    }
  }else{
    supportplot <- NULL
    hlvy_grid_interval <- NULL
  }
  
  result <- list(parameter_space = parameter_space,
                 tree = tree,
                 modfit = modfit,
                 supportplot = supportplot,
                 climblog_matrix = climblog_matrix,
                 brownian_predictors = brownian_predictors,
                 opt.reg = fit$opt.reg,
                 ev.reg = fit$ev.reg,
                 oupar = oupar,
                 hlvy_grid_interval = hlvy_grid_interval,
                 n.par = n.par,
                 V = fit$V,
                 lineages = lineages,
                 fixed.fact = fixed.fact)
  class(result) <- c("slouch", class(result))
  return(result)
}