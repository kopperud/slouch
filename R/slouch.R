#'SLOUCH: Stochastic Linear Ornstein Uhlenbeck Comparative Hypotheses
#'
#'
#'@importFrom Rcpp evalCpp
#'@useDynLib slouch
"_PACKAGE"





#' Title
#'
#' @param phy an object of class 'phylo', must be rooted.
#' @param species a character vector of species tip labels, typically the "species" column in a data frame. This column needs to be an exact match and same order as phy$tip.label
#' @param half_life_values A vector of candidate phylogenetic half-life values to be evaluated in grid search.
#' @param vy_values A vector of candidate stationary variances for the response trait, to be evaluated in grid search.
#' @param response A numeric vector of a trait to be treated as response variable
#' @param me.response Numeric vector of the observational variances of each response trait. E.g if response is a mean trait value, me.response is the within-species SE of the mean.
#' @param fixed.fact factor of regimes on the terminal edges of the tree, in same order as species. If this is used, phy$node.label needs to be filled with the corresponding internal node regimes, in the order of node numbers (root: n+1),(n+2),(n+3), ...
#' @param fixed.cov Direct effect independent variables
#' @param me.fixed.cov Observational variances for direct effect independent variables
#' @param mecov.fixed.cov .
#' @param random.cov Independent variables each modeled as a brownian motion
#' @param me.random.cov Observational variances for the brownian covariates
#' @param mecov.random.cov .
#' @param estimate.Ya Boolean. If true, the intercept K = 1 is expanded to Ya = exp(-a*T.term) and b0 = 1-exp(-a*T.term). If models with categorical covariates are used, this will instead estimate a separate primary optimum for the root niche, "Ya". This only makes sense for non-ultrametric trees. If the tree is ultrametric, the model matrix becomes singular. Defaults to false.
#' @param estimate.bXa Boolean. If true, bXa = bXa = 1-exp(-a*T.term) - (1-(1-exp(-a*T.term))/(a*T.term)) is added to the model matrix, estimating b*Xa. Same requirements as for estimating Ya. Defaults to false.
#' @param support A scalar indicating the size of the support set, defaults to 2 units of log-likelihood.
#' @param convergence Threshold of iterative GLS estimation, when beta is considered converged.
#' @param nCores Use multiple CPU cores in grid-search. If 2 or more cores are used, all print statements are silenced during grid search. Optional, defaults to 1.
#' @param hillclimb TRUE/FALSE whether to use hillclimb parameter estimation routine.
#' @param hillclimb_start Numeric vector of length 2, c(hl, vy), to specify where the hillclimber routine starts. Optional.
#' @param lower lower bounds for the optimization routine, defaults to c(0,0). First entry in vector is half-life, second is stationary variance. When running direct effect models without observational error, it may be useful to specify a positive lower bounds for the stationary variance, e.g c(0, 0.001), since the residual variance-covariance matrix is degenerate when sigma = 0.
#' @param upper upper bounds for the optimization routine, defaults to positive infinite.
#' @param verbose if TRUE, prints each iteration of parameter search. Default FALSE.
#'
#' @return An object of class 'slouch'
#' @export
model.fit.dev2<-function(phy,
                         species = NULL,
                         half_life_values = NULL, 
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
                         upper = Inf,
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
  stopifnot(is.rooted(phy))
  stopifnot((is.numeric(hillclimb_start) & length(hillclimb_start) == 2) | is.null(hillclimb_start))
  if((is.null(half_life_values) | is.null(vy_values)) & !hillclimb){
    stop("Choose at minimum a 1x1 grid, or use the hillclimber routine.")
  }

  # SET DEFAULTS IF NOT SPECIFIED
  if(is.null(me.response)){
    me.response<-diag(rep(0, times=length(response[!is.na(response)])))
  }else{
    me.response<-diag(me.response[!is.na(me.response)])
  }
  
  # SPECIFY COMPONENTS THAT ARE COMMON TO ALL MODELS
  
  Y <- response[!is.na(response)]
  N <- length(Y)
  
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
    lineages <- lapply(1:N, function(e) lineage.constructor(phy, e, regimes)) #; names(lineages) <- phy$tip.label
  }else{
    regimes <- lineages <- NULL
  }

  mrca1 <- mrca(phy)
  times <- node.depth.edgelength(phy)
  ta <- matrix(times[mrca1], nrow=N, dimnames = list(phy$tip.label, phy$tip.label))
  T.term <- times[1:N]
  tia <- times[1:N] - ta
  tja <- t(tia)
  tij <- tja + tia
  
  h.lives<-matrix(data=0, nrow=length(half_life_values), ncol=length(vy_values))
  half_life_values<-rev(half_life_values)
  
  if(!is.null(fixed.cov)){
    if(ncol(as.matrix(fixed.cov))==1) {
      names.fixed.cov <- deparse(substitute(fixed.cov))
    }else{
      stopifnot(!is.null(colnames(fixed.cov)))
      names.fixed.cov <- colnames(fixed.cov)
    }
  }else{
    names.fixed.cov <- NULL
  }
  
  if(!is.null(random.cov)){
    if(ncol(as.matrix(random.cov)) == 1) {
      names.random.cov <- deparse(substitute(random.cov))
    }
    else {
      stopifnot(!is.null(colnames(random.cov)))
      names.random.cov <- colnames(random.cov)
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
               me.fixed.cov = me.fixed.cov,
               mecov.fixed.cov = mecov.fixed.cov,
               random.cov = random.cov,
               me.random.cov = me.random.cov,
               mecov.random.cov = mecov.random.cov,
               Y = Y,
               names.fixed.cov = names.fixed.cov,
               names.random.cov = names.random.cov)
  
  control <- list(verbose = verbose,
                  estimate.Ya = estimate.Ya,
                  estimate.bXa = estimate.bXa,
                  support = support,
                  convergence = convergence)
  
  seed <- ols.seed(tree, pars, control)
  coef.names <- colnames(slouch.modelmatrix(a = 1, hl = 1, tree, pars, control, seed, is.opt.reg = TRUE))
  
  if (is.null(half_life_values)){
    half_life_values <- runif(1, 0, max(times))
  }
  if (is.null(vy_values)){
    vy_values <- runif(1, 0, var(na.exclude(response)))
  }

  if(verbose){
    message("GRID SEARCH PARAMETER SUPPORT")
    cat(c("     hl     ", "vy    ", "support", c(coef.names), "\n"))
  }
  
  #############
  vector_hl_vy <- cbind(sort(rep(half_life_values, length(vy_values)), decreasing = TRUE), rep(vy_values, length(half_life_values)))
  time0 <- Sys.time()
  
  if(nCores > 1 & nCores %% 1 == 0){
    if(!"package:parallel" %in% search()){
      stop("For multicore, please load package with library(parallel)")
    }
    if(.Platform$OS.type == "unix"){
      list_hl_vy <- lapply(seq_len(nrow(vector_hl_vy)), function(e) vector_hl_vy[e,])
      grid_support <- mclapply(list_hl_vy, 
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
  ml_grid <- max(na.exclude(sup2_grid))
  
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
    hl_vy_est <- optim(
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
  
  gof <- matrix(sup2_grid, ncol=length(vy_values), byrow=TRUE, dimnames = list(half_life_values, vy_values))
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
    brownian_predictors <- matrix(data=rbind(seed$theta.X, seed$sigma_squared), nrow=2, ncol=seed$n.pred, dimnames=list(c("Predictor theta", "Predictor variance"), if(seed$n.pred==1) deparse(substitute(random.cov)) else colnames(random.cov)))
  }else{
    brownian_predictors <- NULL
  }
  
  n.par <- length(coef.names) + 2
  modfit<-matrix(data=0, nrow=7, ncol=1, dimnames=list(c("Support", "AIC", "AICc", "SIC", "r squared", "SST", "SSE"),("Value")))
  modfit[1,1] = ml
  modfit[2,1] = -2*ml+2*n.par
  modfit[3,1] = modfit[2,1]+(2*n.par*(n.par+1))/(N-n.par-1)
  modfit[4,1] = -2*ml+log(N)*n.par
  modfit[5,1] = fit$r.squared*100
  modfit[6,1] = fit$sst
  modfit[7,1] = fit$sse
  
  
  if(length(half_life_values) > 1 && length(vy_values) > 1){
    ## All hl + vy in the support interval
    hlsupport <- ifelse(sup2_grid <= ml - support, NA, sapply(grid_support, function(e) e$hl_vy[1]))
    vysupport <- ifelse(sup2_grid <= ml - support, NA, sapply(grid_support, function(e) e$hl_vy[2]))
    
    hlvy_grid_interval <- matrix(c(min(hlsupport, na.rm = TRUE), min(vysupport, na.rm = TRUE),
                                   max(hlsupport, na.rm = TRUE), max(vysupport, na.rm = TRUE)),
                                 ncol = 2, nrow = 2,
                                 dimnames = list(c("Phylogenetic half-life", "Stationary variance"), c("Minimum", "Maximum")))
    
    if(!(all(is.na(gof) | is.infinite(gof)))){
      h.lives <- matrix(0, nrow=length(half_life_values), ncol=length(vy_values), dimnames = list(rev(half_life_values), vy_values))
      for(i in 1:length(vy_values)){
        h.lives[,i]=rev(gof[,i])
      }
      
      supportplot = list(hl = rev(half_life_values),
                         vy = vy_values,
                         z = h.lives)
    }else{
      warning("All support values in grid either NA or +/-Inf - Can't plot.")
    }
  }else{
    supportplot <- NULL
    hlvy_grid_interval <- NULL
  }
  
  print("debug: model.fit.dev2 - slouch in development, use at own risk")
  
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
                 lineages = lineages)
  class(result) <- c("slouch", class(result))
  return(result)
}