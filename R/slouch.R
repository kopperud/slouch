#'SLOUCH: Stochastic Linear Ornstein Uhlenbeck Comparative Hypotheses
#'
#'
"_PACKAGE"



#' Title
#'
#' @param topology
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
#' @param intercept
#' @param ultrametric
#' @param support
#' @param convergence
#' @param plot.angle
#'
#' @return
#' @examples
#' @export
#'

model.fit.dev2<-function(topology, 
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
                        plot.angle = 30,
                        parallel.compute = FALSE,
                        hillclimb = FALSE)
{
  stopifnot(intercept == "root" | is.null(intercept))
  ancestor <- topology
  # SET DEFAULTS IF NOT SPECIFIED
  if(is.null(me.response)){
    me.response<-diag(rep(0, times=length(response[!is.na(response)])))
  }else{
    me.response<-diag(me.response[!is.na(me.response)])
  }

  # SPECIFY COMPONENTS THAT ARE COMMON TO ALL MODELS
  
  Y <- response[!is.na(response)]
  N <- length(Y)
  T.term <- times[terminal.twigs(topology)]
  tia<-tsia(ancestor, time)
  tja<-tsja(ancestor, time)
  term<-terminal.twigs(topology)
  pt<-parse.tree(topology, times)
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
  
  ## Cluster all parameters concerning phylogenetic tree
  treepar <- list(T.term = T.term,
                  tia = tia,
                  tja = tja,
                  term = term,
                  pt = pt,
                  ta = ta,
                  tij = tij,
                  topology = topology,
                  ancestor = ancestor,
                  times = times,
                  species = species)
  
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
  
  all.closures <- regression.closures(treepar, modelpar, seed)
  
  #############
  vector_hl_vy <- cbind(sort(rep(half_life_values, length(vy_values)), decreasing = TRUE), rep(vy_values, length(half_life_values)))
  time0 <- Sys.time()
  
  
  if(hillclimb){
    hl_vy_est <- optim(c(1,1), 
                 all.closures$slouch.regression,
                 gradsearch = TRUE,
                 lower=0, 
                 method="L-BFGS-B")
    
    #print(hl_vy_est)
    #return(hl_vy_est)
    
    best.estimate <- all.closures$slouch.regression(hl_vy_est$par)
    ml <- (-1)*hl_vy_est$value
  }else{
    if(parallel.compute == TRUE){
      n.cores <- detectCores(all.tests = FALSE, logical = TRUE)
      cl <- makeCluster(getOption("cl.cores", n.cores))
      
      estimates <- parApply(cl, vector_hl_vy, 1, all.closures$slouch.regression)
      
      stopCluster(cl)
    }else{
      
      estimates <- apply(vector_hl_vy, 1, all.closures$slouch.regression)
      
    }

    sup2 <- sapply(estimates, function(e) e$support)
    
    gof <- matrix(sup2, ncol=length(vy_values), byrow=TRUE, dimnames = list(half_life_values, vy_values))
    
    ml <- max(na.exclude(gof))
    gof <- ifelse(gof <= ml-support, ml-support, gof) - ml
    
    ## All hl + vy in the support interval
    hlsupport <- ifelse(sup2 <= ml - support, NA, sapply(estimates, function(e) e$hl.est))
    vysupport <- ifelse(sup2 <= ml - support, NA, sapply(estimates, function(e) e$vy.est))
    
    hlvy_grid_interval <- matrix(c(min(hlsupport, na.rm = TRUE), min(vysupport, na.rm = TRUE),
                                   max(hlsupport, na.rm = TRUE), max(vysupport, na.rm = TRUE)),
                                 ncol = 2, nrow = 2,
                                 dimnames = list(c("Phylogenetic half-life", "Stationary variance"), c("Minimum", "Maximum")))
    
    ## Find the regression for which the support value is maximized
    best.estimate <- estimates[[which.max(sup2)]]
  }
  print(paste0("Parameter search done after ",round((Sys.time() - time0), 3)," s."))
  
  
  ############################
  ###### PASTED IN FROM rREG 

  V.est <- best.estimate$V
  V.inverse <- solve(V.est)
  beta1.est <- beta1 <-  best.estimate$beta1
  beta1.var.est <- beta.i.var <- best.estimate$beta1.var
  X <- best.estimate$X
  alpha.est <- best.estimate$alpha.est
  vy.est <- best.estimate$vy.est
  #Y <- best.estimate$Y
  
  ## Calculate evolutionary regression coefficients
  if(!is.null(random.cov) & is.null(fixed.fact)){
    X1 <- calc.X(a = alpha.est, hl = log(2)/alpha.est, treepar, modelpar, seed, is.opt.reg = FALSE)
    ev.beta.i.var<-pseudoinverse(t(X1)%*%V.inverse%*%X1)
    ev.beta.i<-ev.beta.i.var%*%(t(X1)%*%V.inverse%*%Y)
  }else{
    ev.beta.i.var <- NULL
    ev.beta.i <- NULL
  }


  
  opt.reg <- data.frame(cbind(beta1.est, sqrt(diag(beta1.var.est))))
  row.names(opt.reg) <- coef.names
  colnames(opt.reg) <- c("Estimate", "Std. Error")

  message("Model parameters")
  oupar <- matrix(c(alpha.est, log(2)/alpha.est, vy.est), ncol=1, dimnames=list(c("Rate of adaptation", "Phylogenetic half-life", "Stationary variance"), "Estimate"))
  print(oupar)
  
  if(!hillclimb){
    message("Interval of parameters in 3d plot")
    print(hlvy_grid_interval)
  }
  
  message("Optimal regression, estimates + SE")
  print(opt.reg)
  
  if(!is.null(random.cov) & is.null(fixed.fact)){
    ev.reg <- data.frame(cbind(ev.beta.i, sqrt(diag(ev.beta.i.var))))
    row.names(ev.reg) <- coef.names
    colnames(ev.reg) <- c("Estimate", "Std. Error")
    
    message("Ev regr")
    print(ev.reg)
  }
  if(!is.null(random.cov)){
    message("Stochastic predictor")
    print(matrix(data=rbind(seed$theta.X, seed$s.X), nrow=2, ncol=seed$n.pred, dimnames=list(c("Predictor theta", "Predictor variance"), if(seed$n.pred==1) deparse(substitute(random.cov)) else colnames(random.cov))))
  }

  # Model fit statistics
  pred.mean<-X%*%beta1.est
  g.mean<-(t(rep(1, times=N))%*%solve(V.est)%*%Y)/sum(solve(V.est))
  sst<-t(Y-g.mean)%*% solve(V.est)%*%(Y-g.mean)
  sse<-t(Y-pred.mean)%*%solve(V.est)%*%(Y-pred.mean)
  r.squared<-(sst-sse)/sst
  
  n.par<-length(beta1)
  modfit<-matrix(data=0, nrow=7, ncol=1, dimnames=list(c("Support", "AIC", "AICc", "SIC", "r squared", "SST", "SSE"),("Value")))
  modfit[1,1]=ml
  modfit[2,1]=-2*ml+2*(2+n.par)
  modfit[3,1]=modfit[2,1]+(2*(2+n.par)*((2+n.par)+1))/(N-(2+n.par)-1)
  modfit[4,1]=-2*ml+log(N)*(2+n.par)
  modfit[5,1]=r.squared*100
  modfit[6,1]=sst
  modfit[7,1]=sse
  
  print(modfit)
  

  if(!hillclimb){
    # PLOT THE SUPPORT SURFACE FOR HALF-LIVES AND VY
    if(length(half_life_values) > 1 && length(vy_values) > 1){
      h.lives <- matrix(0, nrow=length(half_life_values), ncol=length(vy_values))
      for(i in 1:length(vy_values)){
        h.lives[,i]=rev(gof[,i])
      }
      z<-h.lives
      x<-rev(half_life_values)
      y<-vy_values
      op <- par(bg = "white")
      
      persp(x, y, z, theta = plot.angle, phi = 30, expand = 0.5, col = "NA",
            ltheta = 120, shade = 0.75, ticktype = "detailed",
            xlab = "half-life", ylab = "vy", zlab = "log-likelihood")
    }
  }

  print("debug: model.fit.dev2 w closures")
} # END OF MODEL FITTING FUNCTION