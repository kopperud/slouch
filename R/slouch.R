#'SLOUCH: Stochastic Linear Ornstein Uhlenbeck Comparative Hypotheses
#'
#'
"_PACKAGE"





#' Title
#'
#' @param ancestor Vector of ancestors for each given node in tree
#' @param times A numeric vector; time from root to given node
#' @param half_life_values A vector of candidate phylogenetic half-life values to be evaluated in grid search.
#' @param vy_values A vector of candidate stationary variances for the response trait, to be evaluated in grid search.
#' @param response A numeric vector of a trait to be treated as response variable
#' @param me.response Numeric vector of the observational variances of each response trait. E.g if response is a mean trait value, me.response is the within-species SE of the mean.
#' @param fixed.fact Factor
#' @param fixed.cov Direct effect independent variables
#' @param me.fixed.cov Observational variances for direct effect independent variables
#' @param mecov.fixed.cov .
#' @param random.cov Independent variables each modeled as a brownian motion
#' @param me.random.cov Observational variances for the brownian covariates
#' @param mecov.random.cov .
#' @param ultrametric Deprecated
#' @param intercept NULL or "root". If NULL, model matrix is expanded with Ya for categorical models, or b0 + bXa for models without categorical variables.
#' @param support A scalar indicating the size of the support set, defaults to 2 units of log-likelihood.
#' @param convergence Threshold of iterative GLS estimation, when beta is considered converged.
#' @param multicore Use of multiple CPU cores. Using multicore, prints are silenced during parameter search.
#' @param ncores Number of CPU cores to be used, optional.
#' @param hillclimb TRUE/FALSE whether to use hillclimb parameter estimation routine.
#' @param hillclimb_start Numeric vector of length 2, c(hl, vy), to specify where the hillclimber routine starts. Optional.
#' @param verbose if TRUE, prints each iteration of parameter search. Default FALSE.
#'
#' @return An object of class 'slouch'
#' @export
model.fit.dev2<-function(ancestor, 
                         times, 
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
                         ultrametric=TRUE, 
                         intercept= if(ultrametric==TRUE) "root" else NULL, 
                         support = 2, 
                         convergence = 0.000001,
                         multicore = FALSE,
                         ncores = NULL,
                         hillclimb = FALSE,
                         hillclimb_start = NULL,
                         verbose = FALSE)
{
  stopifnot(intercept == "root" | is.null(intercept))
  stopifnot((is.numeric(hillclimb_start) & length(hillclimb_start) == 2) | is.null(hillclimb_start))
  ancestor <- ancestor
  # SET DEFAULTS IF NOT SPECIFIED
  if(is.null(me.response)){
    me.response<-diag(rep(0, times=length(response[!is.na(response)])))
  }else{
    me.response<-diag(me.response[!is.na(me.response)])
  }
  
  # SPECIFY COMPONENTS THAT ARE COMMON TO ALL MODELS
  
  Y <- response[!is.na(response)]
  N <- length(Y)
  T.term <- times[terminal.twigs(ancestor)]
  tia<-tsia(ancestor, times)
  tja<-tsja(ancestor, times)
  
  term<-terminal.twigs(ancestor)
  pt<-parse.tree(ancestor, times)
  ta<-pt$bt
  tij<-pt$dm
  
  h.lives<-matrix(data=0, nrow=length(half_life_values), ncol=length(vy_values))
  half_life_values<-rev(half_life_values)
  
  ## Establish beta names and descriptor
  #factor.names <- coef.names.factor(fixed.fact, random.cov, fixed.cov, intercept)
  
  names.fixed.cov <- if(!is.null(fixed.cov)){
    if(ncol(as.matrix(fixed.cov))==1) deparse(substitute(fixed.cov)) else colnames(fixed.cov)
  }else{
    NULL
  }
  
  names.random.cov <- if(!is.null(random.cov)){
    if(ncol(as.matrix(random.cov)) == 1) deparse(substitute(random.cov)) else colnames(random.cov)
  }else{
    NULL
  }
  
  if(!is.null(fixed.fact)){
    regime.specs <- as.factor(fixed.fact)
    regimes1 <- regimes(ancestor, times, regime.specs, term)
    epochs1 <- epochs(ancestor, times, term)
  }else{
    regimes1 <- NULL
    epochs1 <- NULL
  }
  
  ## Cluster all parameters concerning phylogenetic tree
  treepar <- list(T.term = T.term,
                  tia = tia,
                  tja = tja,
                  term = term,
                  pt = pt,
                  ta = ta,
                  tij = tij,
                  ancestor = ancestor,
                  times = times,
                  #species = species,
                  regimes1 = regimes1,
                  epochs1 = epochs1)
  
  ## Cluster parameters concerning the type of model being run
  modelpar <- list(response = response,
                   me.response = me.response,
                   fixed.fact = fixed.fact,
                   fixed.cov = fixed.cov,
                   me.fixed.cov = me.fixed.cov,
                   mecov.fixed.cov = mecov.fixed.cov,
                   random.cov = random.cov,
                   me.random.cov = me.random.cov,
                   mecov.random.cov = mecov.random.cov,
                   Y = Y,
                   N = N,
                   h.lives = h.lives,
                   half_life_values = half_life_values,
                   vy_values = vy_values,
                   support = support,
                   convergence = convergence,
                   intercept = intercept,
                   regime.specs = as.factor(fixed.fact),
                   factor.exists = !is.null(fixed.fact),
                   names.fixed.cov = names.fixed.cov,
                   names.random.cov = names.random.cov,
                   verbose = verbose)
  
  
  
  seed <- ols.seed(treepar, modelpar)
  coef.names <- colnames(calc.X(a = 1, hl = 1, treepar, modelpar, seed, is.opt.reg = TRUE))
  
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
  
  if(multicore == TRUE){
    if(!"package:parallel" %in% search()){
      stop("For multicore, please load package with library(parallel)")
    }
    if(is.null(ncores)){
      ncores <- detectCores(all.tests = FALSE, logical = TRUE)
    }
    if(.Platform$OS.type == "unix"){
      stop("Parallel computing parameter search is not supported on unix-based systems. To be fixed.")
      list_hl_vy <- unname(split(vector_hl_vy, rep(1:nrow(vector_hl_vy), each = ncol(vector_hl_vy))))
      grid_support <- mclapply(list_hl_vy, reg, modelpar, treepar, seed, mc.cleanup = TRUE, mc.cores = ncores)
    }else{
      cl <- parallel::makeCluster(getOption("cl.cores", ncores))
      parallel::setDefaultCluster(cl)
      parallel::clusterExport(cl, c("modelpar", "treepar", "seed"), envir = environment())
      grid_support <- parallel::parApply(cl, vector_hl_vy, 1, function(e) reg(e, modelpar, treepar, seed))
      parallel::stopCluster(cl)
    }
  }else{
    grid_support <- apply(vector_hl_vy, 1, reg, modelpar, treepar, seed)
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
      fn = function(e, ...){hcenv$k <- hcenv$k +1; tmp <- reg(e, modelpar, treepar, seed, ...); hcenv$climblog[[toString(hcenv$k)]] <- tmp; return((-1)*tmp$support) },
      gridsearch = TRUE,
      lower=0, 
      method="L-BFGS-B")
    
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
  fit <- reg(besthl_vy, modelpar, treepar, seed, gridsearch=FALSE)
  if(verbose){
    print(paste0("Parameter search done after ",round((Sys.time() - time0), 3)," s."))
  }
  
  ############################
  
  alpha <- log(2) / fit$hl_vy[1]
  
  oupar <- matrix(c(log(2) / fit$hl_vy[1], 
                    fit$hl_vy[1], 
                    fit$hl_vy[2], 
                    mean((1-(1-exp(-alpha*T.term))/(alpha*T.term)))), 
                  ncol=1, 
                  dimnames=list(c("Rate of adaptation", "Phylogenetic half-life", "Stationary variance", "Phylogenetic correction factor"), "Estimate"))
  
  if(!is.null(random.cov)){
    brownian_predictors <- matrix(data=rbind(seed$theta.X, seed$s.X), nrow=2, ncol=seed$n.pred, dimnames=list(c("Predictor theta", "Predictor variance"), if(seed$n.pred==1) deparse(substitute(random.cov)) else colnames(random.cov)))
  }else{
    brownian_predictors <- NULL
  }
  
  n.par <- length(fit$opt.reg$coefficients[,1]) + 2
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
      # PLOT THE SUPPORT SURFACE FOR HALF-LIVES AND VY
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
                 tree = list(tia = tia,
                             tja = tja,
                             tij = tij,
                             ta = ta),
                 modfit = modfit,
                 supportplot = supportplot,
                 climblog_matrix = climblog_matrix,
                 brownian_predictors = brownian_predictors,
                 opt.reg = fit$opt.reg,
                 ev.reg = fit$ev.reg,
                 oupar = oupar,
                 hlvy_grid_interval = hlvy_grid_interval,
                 n.par = n.par,
                 V = fit$V)
  class(result) <- c("slouch", class(result))
  return(result)
}