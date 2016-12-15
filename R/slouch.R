#'SLOUCH: Stochastic Linear Ornstein Uhlenbeck Comparative Hypotheses
#'
#'
"_PACKAGE"




#' Title
#'
#' @param ancestor 
#' @param times 
#' @param half_life_values 
#' @param vy_values 
#' @param response 
#' @param me.response 
#' @param fixed.fact 
#' @param fixed.cov 
#' @param me.fixed.cov 
#' @param mecov.fixed.cov 
#' @param random.cov 
#' @param me.random.cov 
#' @param mecov.random.cov 
#' @param ultrametric 
#' @param intercept 
#' @param support 
#' @param convergence 
#' @param multicore 
#' @param ncores 
#' @param hillclimb 
#' @param hillclimb_start 
#'
#' @return
#' @export
#'
#' @examples
model.fit.dev2<-function(ancestor, 
                        times, 
                        half_life_values, 
                        vy_values, 
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
                        hillclimb_start = c(runif(1, 0, max(times)), runif(1, 0, var(na.exclude(response))))) # runif(2,0,1)
{
  stopifnot(intercept == "root" | is.null(intercept))
  stopifnot(is.numeric(hillclimb_start) & length(hillclimb_start) == 2)
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
                  #ancestor = ancestor,
                  ancestor = ancestor,
                  times = times,
                  species = species,
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
                   names.random.cov = names.random.cov)
  

  
  seed <- ols.seed(treepar, modelpar)
  coef.names <- colnames(calc.X(a = 1, hl = 1, treepar, modelpar, seed, is.opt.reg = TRUE))

  

  message("GRID SEARCH PARAMETER SUPPORT")
  cat(c("     hl     ", "vy    ", "support", c(coef.names), "\n"))
  
  #############
  vector_hl_vy <- cbind(sort(rep(half_life_values, length(vy_values)), decreasing = TRUE), rep(vy_values, length(half_life_values)))
  time0 <- Sys.time()
  
  
  if(hillclimb){
    hcenv <- environment()
    hcenv$k <- 0
    climblog <- list()
    hl_vy_est <- optim(#c(1,1), 
                 hillclimb_start,
                 function(e, ...){hcenv$k <- hcenv$k +1; tmp <- reg(e, modelpar, treepar, seed, ...); hcenv$climblog[[toString(hcenv$k)]] <- tmp; return((-1)*tmp$support) },
                 gridsearch = TRUE,
                 lower=0, 
                 method="L-BFGS-B")
    besthl_vy <- hl_vy_est$par
    ml <- (-1)*hl_vy_est$value
    
    ## Matrix for plotting the route of hillclimber
    climblog_matrix <- data.frame(index = 1:length(climblog), hl = sapply(climblog, function(e) e$hl_vy[[1]]), vy = sapply(climblog, function(e) e$hl_vy[2]), loglik = sapply(climblog, function(e) e$support))
    hlvy_grid_interval <- NULL
    
  }else{
    climblog_matrix <- NULL
    
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
        #print(grid_support)
      }else{
        cl <- parallel::makeCluster(getOption("cl.cores", ncores))
        parallel::setDefaultCluster(cl)
        #grid_support <- parallel::parApply(cl, vector_hl_vy, 1, function(e) reg(e, modelpar, treepar, seed))
        parallel::clusterExport(cl, c("modelpar", "treepar", "seed"), envir = environment())
        grid_support <- parallel::parApply(cl, vector_hl_vy, 1, function(e) reg(e, modelpar, treepar, seed))
        parallel::stopCluster(cl)
      }
    }else{
      
      grid_support <- apply(vector_hl_vy, 1, reg, modelpar, treepar, seed)
      
    }

    sup2 <- sapply(grid_support, function(e) e$support)
    ml <- max(na.exclude(sup2))
    
    gof <- matrix(sup2, ncol=length(vy_values), byrow=TRUE, dimnames = list(half_life_values, vy_values))
    
    gof <- ifelse(gof <= ml-support, ml-support, gof) - ml
    
    ## All hl + vy in the support interval
    hlsupport <- ifelse(sup2 <= ml - support, NA, sapply(grid_support, function(e) e$hl_vy[1]))
    vysupport <- ifelse(sup2 <= ml - support, NA, sapply(grid_support, function(e) e$hl_vy[2]))
    
    hlvy_grid_interval <- matrix(c(min(hlsupport, na.rm = TRUE), min(vysupport, na.rm = TRUE),
                                   max(hlsupport, na.rm = TRUE), max(vysupport, na.rm = TRUE)),
                                 ncol = 2, nrow = 2,
                                 dimnames = list(c("Phylogenetic half-life", "Stationary variance"), c("Minimum", "Maximum")))
    
    ## Find the regression for which the support value is maximized
    besthl_vy = vector_hl_vy[which.max(sup2),]
  }
  best.estimate <- reg(besthl_vy, modelpar, treepar, seed, gridsearch=FALSE)
  
  print(paste0("Parameter search done after ",round((Sys.time() - time0), 3)," s."))
  

  
  ############################
  ###### PASTED IN FROM rREG 
  V <- best.estimate$V
  ev.reg <- best.estimate$ev.reg
  opt.reg <- best.estimate$opt.reg
  r.squared <- best.estimate$r.squared
  sse <- best.estimate$sse
  sst <- best.estimate$sst
  alpha.est <- log(2) / best.estimate$hl_vy[1]
  vy.est <- best.estimate$hl_vy[2]

  oupar <- matrix(c(alpha.est, 
                    log(2)/alpha.est, 
                    vy.est, 
                    mean((1-(1-exp(-alpha.est*T.term))/(alpha.est*T.term)))), 
                  ncol=1, 
                  dimnames=list(c("Rate of adaptation", "Phylogenetic half-life", "Stationary variance", "Phylogenetic correction factor"), "Estimate"))

  if(!is.null(random.cov)){
    brownian_predictors <- matrix(data=rbind(seed$theta.X, seed$s.X), nrow=2, ncol=seed$n.pred, dimnames=list(c("Predictor theta", "Predictor variance"), if(seed$n.pred==1) deparse(substitute(random.cov)) else colnames(random.cov)))
  }else{
    brownian_predictors <- NULL
  }

  n.par <- length(opt.reg$coefficients[,1]) + 2
  modfit<-matrix(data=0, nrow=7, ncol=1, dimnames=list(c("Support", "AIC", "AICc", "SIC", "r squared", "SST", "SSE"),("Value")))
  modfit[1,1]=ml
  modfit[2,1]=-2*ml+2*n.par
  modfit[3,1]=modfit[2,1]+(2*n.par*(n.par+1))/(N-n.par-1)
  modfit[4,1]=-2*ml+log(N)*n.par
  modfit[5,1]=r.squared*100
  modfit[6,1]=sst
  modfit[7,1]=sse

  supportplot <- NULL
  if(!hillclimb){
    if(length(half_life_values) > 1 && length(vy_values) > 1){
      if(!(all(is.na(gof) | is.infinite(gof)))){
        
        # PLOT THE SUPPORT SURFACE FOR HALF-LIVES AND VY
        h.lives <- matrix(0, nrow=length(half_life_values), ncol=length(vy_values))
        for(i in 1:length(vy_values)){
          h.lives[,i]=rev(gof[,i])
        }
        z<-h.lives
        x<-rev(half_life_values)
        y<-vy_values
        op <- par(bg = "white")
        
        # persp(x, y, z, theta = plot.angle, phi = 30, expand = 0.5, col = "NA",
        #       ltheta = 120, shade = 0.75, ticktype = "detailed",
        #       xlab = "half-life", ylab = "vy", zlab = "log-likelihood")
        supportplot = list(x = x,
                           y = y,
                           z = z)
      }else{
        warning("All support values in grid either NA or +/-Inf - Can't plot.")
      }
    }
  }
  
  print("debug: model.fit.dev2 - slouch in development, use at own risk")
  
  result <- list(grid_support =  if (!hillclimb) grid_support else climblog,
                 tree = list(tia = tia,
                             tja = tja,
                             tij = tij,
                             ta = ta),
                 modfit = modfit,
                 supportplot = supportplot,
                 climblog_matrix = climblog_matrix,
                 brownian_predictors = brownian_predictors,
                 opt.reg = opt.reg,
                 ev.reg = ev.reg,
                 oupar = oupar,
                 hlvy_grid_interval = hlvy_grid_interval,
                 n.par = n.par,
                 V = V)
  class(result) <- c("slouch", class(result))
  return(result)
}