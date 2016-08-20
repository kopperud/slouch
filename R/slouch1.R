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
                        intercept="root", 
                        ultrametric=TRUE, 
                        support=NULL, 
                        convergence=NULL, 
                        plot.angle=30, 
                        parallel.compute = FALSE)
{
  ancestor <- topology
  # SET DEFAULTS IF NOT SPECIFIED
  
  if(is.null(support)) support <- 2
  if(is.null(convergence)) convergence <- 0.000001
  if(is.null(me.response)){
    me.response<-diag(rep(0, times=length(response[!is.na(response)])))
  }else{
    me.response<-diag(me.response[!is.na(me.response)])
  }
  
  
  # DETERMINE MODEL STRUCTURE FROM INPUT AND WRITE A SUMMARY TO THE R CONSOLE
  
  if(is.null(fixed.fact) && is.null(fixed.cov) && is.null(random.cov)) model.type <- "IntcptReg";
  if(!is.null(fixed.fact) && is.null(fixed.cov) && is.null(random.cov)) model.type <- "ffANOVA";
  if(!is.null(fixed.fact) && !is.null(fixed.cov) && is.null(random.cov)) model.type <-"ffANCOVA";
  if(!is.null(fixed.fact) && is.null(fixed.cov) && !is.null(random.cov)) model.type <- "mmANCOVA";
  if(!is.null(fixed.fact) && !is.null(fixed.cov) && !is.null(random.cov)) model.type <- "mmfANCOVA";
  if(is.null(fixed.fact) && is.null(fixed.cov) && !is.null(random.cov)) model.type <- "rReg";
  if(is.null(fixed.fact) && !is.null(fixed.cov) && is.null(random.cov)) model.type <- "fReg";
  if(is.null(fixed.fact) && !is.null(fixed.cov) && !is.null(random.cov)) model.type <- "mfReg";
  
  # Write type of model to screen
  message("")
  message("MODEL SUMMARY")
  message("")
  if(model.type=="IntcptReg")
  {
    message("You have specified an OU model for a response variable regressed on a grand mean, i.e. one global optima");
    if(ultrametric==FALSE)
    {
      GS_head<-c("Ya", "Theta_Global")
      #n.par<-2 # Bjorn: unused?
    }
    else
    {
      GS_head<-("Theta_Global")
      #n.par<-1 # Bjorn: unused?
    }
  }
  else
    if(model.type=="ffANOVA" )
    {
      message("You have specified an OU model for a response variable modeled on optima determined by fixed, categorical predictor variables");
      if(is.null(intercept)) GS_head<-c("Ya", levels(as.factor(fixed.fact))) else GS_head<-levels(as.factor(fixed.fact));
    }
  
  else
    if(model.type=="ffANCOVA")
    {
      message("You have specified an OU model for a response variable modeled on optima determined by both fixed categorical predictors and an instantaneous scaling with a fixed covariate");
      if(is.null(intercept)) GS_head<-c("Ya", levels(as.factor(fixed.fact))) else GS_head<-levels(as.factor(fixed.fact));
      
    }
  
  
  else
    
    if(model.type=="mmANCOVA")
    {
      message("You have specified an OU model for a response variable modeled on optima determined by both fixed, categorical factors as well as covariates which themselves randomly evolve (modeled as Brownian-motions)");
      
      
      if(is.null(intercept)) GS_head<-c("Ya", levels(as.factor(fixed.fact))) else GS_head<-levels(as.factor(fixed.fact));
    }
  
  if(model.type=="mmfANCOVA")
  {
    message("You have specified an OU model for a response variable modeled on optima determined by both fixed, categorical factors as well as covariates which themselves randomly evolve (modeled as Brownian-motions)");
    
    
    if(is.null(intercept)) GS_head<-c("Ya", levels(as.factor(fixed.fact))) else GS_head<-levels(as.factor(fixed.fact));
  }
  else
    
    
    
    if(model.type=="rReg") message("You have specified an OU model for a response variable modeled on optima that are determined by randomly evolving covariates (modeled as Brownian-motions)")
  
  else
    if(model.type=="fReg") message("You have specified an OU model for a response variable modeled on optima that are determined by an instantaneous scaling with fixed covariates")
  
  else
    if(model.type=="mfReg") message("You have specified an OU model for a response variable modeled on optima that are determined by both an instantaneous scaling with fixed covariates and randomly evolving covariates (modeled as Brownian-motions)");
  message("")
  
  # Summarize dataset, response, predictors,  tree height and sample size and write to screen
  
  ms<-list(Dataset=search()[2], Response=deparse(substitute(response)), Fixed.factor=deparse(substitute(fixed.fact)),Fixed.covariates=deparse(substitute(fixed.cov)), Random.covariates=deparse(substitute(random.cov)), Sample.size=length(response[!is.na(response)]), Tree.height=max(times), Model.type=substitute(model.type))
  ms<-as.matrix(ms)
  colnames(ms)<-"Summary"
  print(ms)
  message("")
  message("GRID SEARCH PARAMETER SUPPORT")
  message("")
  
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
  
  ## Cluster all parameters concerning phylogenetic tree
  treepar <- list(T.term = T.term,
                  tia = tia,
                  tja = tja,
                  term = term,
                  pt = pt,
                  ta = ta,
                  tij = tij,
                  ultrametric = ultrametric,
                  topology = topology,
                  ancestor = ancestor,
                  times = times,
                  species = species)
  
  ## Cluster parameters concerning the type of model being run
  modelpar <- list(model.type = model.type,
                   response = response,
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
                   factor.exists = !is.null(fixed.fact))
  
  ## Create seed with OLS
  #######################
  ## --------------------
  #######################
  
  seed <- ols.seed(treepar, modelpar)
  
  #return(seed)
  
  #######################
  ## --------------------
  #######################
  
  
  
  all.closures <- regression.closures(treepar, modelpar, seed)
  
  
  #list2env(all.closures, envir = environment())
  
  ## Test hl, vy == 1,1
  #slouch.regression(c(1,1))

  vector_hl_vy <- cbind(sort(rep(half_life_values, length(vy_values)), decreasing = TRUE), rep(vy_values, length(half_life_values)))
  estimates <- apply(vector_hl_vy, 1, all.closures$slouch.regression)
  
  sup2 <- sapply(estimates, function(e) e$support)
  gof <- matrix(sup2, ncol=length(vy_values), byrow=TRUE, dimnames = list(half_life_values, vy_values))
  
  
  ml<-max(na.exclude(gof))
  gof <- ifelse(gof <= ml-support, ml-support, gof) - ml
  
  ############################
  ###### PASTED IN FROM rREG 
  
  ## Very bad: don't add list to env
  list2env(seed, envir = environment())
  
  ## Find the regression for which the support value is maximized
  best.estimate <- estimates[[which.max(sup2)]]
  V.est <- best.estimate$V; V.inverse <- solve(V.est)
  beta1.est <- beta1 <-  best.estimate$beta1
  beta1.var.est <- beta.i.var <- best.estimate$beta1.var
  X <- best.estimate$X
  alpha.est <- best.estimate$alpha.est
  vy.est <- best.estimate$vy.est
  Y <- best.estimate$Y
  
  X1<-cbind(1, pred)
  ev.beta.i.var<-pseudoinverse(t(X1)%*%V.inverse%*%X1)
  ev.beta.i<-ev.beta.i.var%*%(t(X1)%*%V.inverse%*%Y)
  
  ## Calculate model fit stats
  pred.mean<-X%*%beta1.est
  g.mean<-(t(rep(1, times=N))%*%solve(V.est)%*%Y)/sum(solve(V.est))
  sst<-t(Y-g.mean)%*% solve(V.est)%*%(Y-g.mean)
  sse<-t(Y-pred.mean)%*%solve(V.est)%*%(Y-pred.mean)
  r.squared<-(sst-sse)/sst
  
  ## Calculate AIC, AICc
  n.par<-length(beta1)
  aic <- -2*ml+2*(2+n.par)
  aicc <- aic +(2*(2+n.par)*((2+n.par)+1))/(N-(2+n.par)-1)

  #print(V.est)
  message("Optimal regression, estimates + SE")
  print(cbind(beta1.est, sqrt(diag(beta1.var.est))))
  message("Ev regr")
  print(cbind(ev.beta.i, sqrt(diag(ev.beta.i.var))))
  message("Alpha, hl, vy")
  print(c(alpha.est, 
          if(alpha.est == Inf) 0 else log(2)/alpha.est, 
          vy.est))
  
  
  modfit<-matrix(data=0, nrow=7, ncol=1, dimnames=list(c("Support", "AIC", "AICc", "SIC", "r squared", "SST", "SSE"),("Value")))
  modfit[1,1]=ml
  modfit[2,1]=-2*ml+2*(2+n.par)
  modfit[3,1]=modfit[2,1]+(2*(2+n.par)*((2+n.par)+1))/(N-(2+n.par)-1)
  modfit[4,1]=-2*ml+log(N)*(2+n.par)
  modfit[5,1]=r.squared*100
  modfit[6,1]=sst
  modfit[7,1]=sse
  
  print(modfit)
  
  # PLOT THE SUPPORT SURFACE FOR HALF-LIVES AND VY
  
  
  if(length(half_life_values) > 1 && length(vy_values) > 1){
      z1<-gof
      for(i in 1:length(vy_values)){
        h.lives[,i]=rev(z1[,i])
      }
      z<-h.lives
      x<-rev(half_life_values)
      y<-vy_values
      op <- par(bg = "white")
      
      persp(x, y, z, theta = plot.angle, phi = 30, expand = 0.5, col = "NA",
            ltheta = 120, shade = 0.75, ticktype = "detailed",
            xlab = "half-life", ylab = "vy", zlab = "log-likelihood")
  }

  print("debug: model.fit.dev w closures")
  
  
} # END OF MODEL FITTING FUNCTION


