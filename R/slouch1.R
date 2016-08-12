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
  
  
  all.closures <- regression.closures(treepar, modelpar, seed)
  
  

  
  
  # PLOT THE SUPPORT SURFACE FOR HALF-LIVES AND VY
  
  
  if(length(half_life_values) > 1 && length(vy_values) > 1){
    #if(length(gof[gof>-2])){
      z1<-gof
      for(i in 1:length(vy_values)){
        h.lives[,i]=rev(z1[,i])
      }
      z<-h.lives
      x<-rev(half_life_values)
      y<-vy_values
      op <- par(bg = "white")
      
      # plot.slouch.x <<- x
      # plot.slouch.y <<- y
      # plot.slouch.loglik <<- z
      
      #persp(x, y, z, theta = plot.angle, phi = 30, expand = 0.5, col = "NA") ## plot.angle = 30 default
      persp(x, y, z, theta = plot.angle, phi = 30, expand = 0.5, col = "NA",
            ltheta = 120, shade = 0.75, ticktype = "detailed",
            xlab = "half-life", ylab = "vy", zlab = "log-likelihood")
    #}else{
      #print(sort(gof, decreasing = TRUE))
    #}
  }
  #plot.coord <- cbind(rep(x, length(y)), rep(rev(y), length(x)))
  
  
  
  
  
  # fixed.cov.names <- if(is.null(dim(fixed.cov))) deparse(substitute(fixed.cov)) else colnames(fixed.cov)
  # random.cov.names <- if(is.null(dim(random.cov))) deparse(substitute(random.cov)) else colnames(random.cov)
  #
  # output5 <- list(model.type = model.type,
  #                 alpha.est = alpha.est,
  #                 vy.est = vy.est,
  #                 pred.mean = pred.mean,
  #                 g.mean = g.mean,
  #                 sst = sst,
  #                 sse = sse,
  #                 r.squared = r.squared,
  #                 ml = ml,
  #                 n.pred = n.pred,
  #                 ultrametric = ultrametric,
  #                 intercept = intercept,
  #                 N = N,
  #                 beta.est = beta.est,
  #                 theta.X = theta.X,
  #                 s.X = s.X,
  #                 beta.var.est = beta.var.est,
  #                 X = X,
  #                 x.ols = x.ols,
  #                 random.cov = random.cov,
  #                 corrected_betas = corrected_betas,
  #                 fixed.fact = fixed.fact,
  #                 random.cov.names = random.cov.names,
  #                 glsyx.beta1 = glsyx.beta1,
  #                 glsyx.beta1.var = glsyx.beta1.var,
  #                 fixed.cov = fixed.cov,
  #                 me.fixed.cov = me.fixed.cov,
  #                 fixed.cov.names = fixed.cov.names)
  #
  # result <- new("slouch", output5)
  # return(result)
  
  
  
  # MODEL OUTPUT
  # alpha, half-lives, correction factor, v
  
  
  message("==================================================")
  half.life<-log(2)/alpha.est
  c.factor<-mean(1-(1-exp(-alpha.est*T))/(alpha.est*T))
  modeloutput<-matrix(data=0, nrow=4, ncol=1, dimnames=list(c("Rate of adaptation ", "Phylogenetic half-life ","Phylogenetic correction factor", "Stationary variance "), "    Estimate"))
  modeloutput[1, 1]=alpha.est; modeloutput[2, 1]=half.life; modeloutput[3,1]=c.factor; modeloutput[4,1]=vy.est;   ##### Rememeber to output s.X
  
  
  
  
  
  modfit<-matrix(data=0, nrow=7, ncol=1, dimnames=list(c("Support", "AIC", "AICc", "SIC", "r squared", "SST", "SSE"),("Value")))
  
  
  #if(ultrametric==TRUE) n.par=1+n.pred else n.par=3+n.pred
  
  if(model.type=="ffANOVA" || model.type=="fReg" || model.type=="ffANCOVA") n.par<-length(gls.beta0)
  if(model.type == "mmANCOVA" || model.type=="rReg" || model.type=="mfReg" || model.type=="mmfANCOVA")   n.par<-length(beta1)
  
  modfit[1,1]=ml
  modfit[2,1]=-2*ml+2*(2+n.par)
  modfit[3,1]=modfit[2,1]+(2*(2+n.par)*((2+n.par)+1))/(N-(2+n.par)-1)
  modfit[4,1]=-2*ml+log(N)*(2+n.par)
  modfit[5,1]=r.squared*100
  modfit[6,1]=sst
  modfit[7,1]=sse
  
  message("");
  message("BEST ESTIMATES & MODEL FIT");message("");
  message("==================================================");
  message("MODEL PARAMETERS");
  print(modeloutput);message("");
  
  
  # predictor means and variances for random predictors
  
  if(model.type == "mmANCOVA" || model.type=="rReg" || model.type=="mfReg" || model.type=="mmfANCOVA")
  {
    print(matrix(data=rbind(theta.X, s.X), nrow=2, ncol=n.pred, dimnames=list(c("Predictor theta", "Predictor variance"), if(n.pred==1) deparse(substitute(random.cov)) else colnames(random.cov))));
    message("");
  }
  
  # PRIMARY OPTIMA OR REGRESSION SLOPE ESTIMATES
  
  message("--------------------------------------------------");
  message("PRIMARY OPTIMA");message("");
  
  
  if(model.type=="IntcptReg")
  {
    if(ultrametric==TRUE || alpha.est==Inf || alpha.est>=1000000000000000){
      Intercept<-matrix(nrow=1, ncol=2, dimnames=list(("Theta_global"), c("Estimate", "Std.error")))
      Intercept[,1]<-gls.beta0
      Intercept[,2]<-sqrt(beta.i.var)}
    else {
      Intercept<-matrix(data=0, nrow=2, ncol=1, dimnames=list(c("Bo", "Ya"), ("     Estimate")))
      Intercept[1,1]<-beta.i[1]
      Intercept[2,1]<-beta.i[2]
    }
    print(Intercept); message("")
  }
  
  if(model.type=="ffANOVA")
  {
    std<-sqrt(diag(beta.i.var))
    
    optima<-matrix(data=0, nrow=ncol(X), ncol=2, dimnames = list(colnames(X), c("Estimates", "Std.error")));
    optima[,1] = gls.beta0;
    optima[,2] = std;
    
    
    reg <- set.of.regimes(topology,regime.specs);
    root.reg<-as.character(regime.specs[times==0])
    nonroot.reg<-as.character(reg[reg != root.reg])
    
    
    if(is.null(intercept))
    {
      if(ncol(X) == length(reg)) message ("The ancestral state (Ya) parameter was dropped from this model as there is not enough information to estimate it")  else
        if(ncol(X)<length(reg)) message ("Ya and the parameter at the root were dropped") else
          message("this model does not drop Ya as it may influence the other parameters")
    }
    else
    {
      if(intercept=="root") message(root.reg, " ", "mapped to the root of the tree and includes the coefficent for the ancestral state (Ya)") else
        message("you set the intercept coefficent to a value of", " ", intercept,". Ya is not the true ancestral state anymore")
    }
    print(optima);message("");
  }
  
  
  if(model.type== "fReg")
  {
    std<-sqrt(diag(beta.i.var))
    
    optima<-matrix(data=0, nrow=(nrow(gls.beta0)), ncol=2, dimnames=list(c("Bo", if(is.null(dim(fixed.cov))) deparse(substitute(fixed.cov)) else colnames(fixed.cov)), c("Estimate", "Std. Error")))
    optima[,1] = gls.beta0;
    optima[,2] = std;
    
    
    corrected_beta_values<-matrix(data= c(corrected_betas), nrow=(nrow(beta1)), ncol=1, dimnames=list(c("Bo", if(is.null(dim(fixed.cov))) deparse(substitute(fixed.cov)) else colnames(fixed.cov)), c("Bias-corr. regression parameters")))
    
    
    print(optima);message("");message("");
    print(corrected_beta_values);message("");
    
  }
  
  if(model.type=="ffANCOVA")
  {
    std<-sqrt(diag(beta.i.var))
    
    
    optima<-matrix(data=0, nrow=ncol(X), ncol=2, dimnames = list(c(as.character(levels(fixed.fact)), if(is.null(dim(fixed.cov))) deparse(substitute(fixed.cov)) else colnames(fixed.cov)), c("Estimates", "Std.error")));
    optima[,1] = gls.beta0;
    optima[,2] = std;
    
    corrected_beta_values<-matrix(data= c(corrected_betas), nrow=(nrow(beta1)), ncol=1, dimnames = list(c(as.character(levels(fixed.fact)), if(is.null(dim(fixed.cov))) deparse(substitute(fixed.cov)) else colnames(fixed.cov)), c("Bias-corr. regression parameters")));
    
    print(optima);message("");message("");
    print(corrected_beta_values);message("");
    
  }
  
  
  
  if(model.type  == "mmANCOVA")
  {
    std<-sqrt(diag(beta.i.var))
    if(length(X[1,]) > length(x.ols[1,])) optima<-matrix(data=0, nrow=ncol(X), ncol=2, dimnames = list(c(c("Ya",as.character(levels(fixed.fact))), if(is.null(dim(random.cov))) deparse(substitute(random.cov)) else colnames(random.cov)), c("Estimates", "Std.error")))
    else
      
      optima<-matrix(data=0, nrow=ncol(X), ncol=2, dimnames = list(c(as.character(levels(fixed.fact)), if(is.null(dim(random.cov))) deparse(substitute(random.cov)) else colnames(random.cov)), c("Estimates", "Std.error")));
    
    corrected_beta_values<-matrix(data= c(corrected_betas), nrow=(nrow(beta1)), ncol=1, dimnames = list(c(as.character(levels(fixed.fact)), if(is.null(dim(random.cov))) deparse(substitute(random.cov)) else colnames(random.cov)), c("Bias-corr. regression parameters")));
    
    
    optima[,1] = gls.beta1;
    optima[,2] = std;
    print(optima);message("");message("");
    print(corrected_beta_values);message("");
  }
  
  if(model.type  == "mmfANCOVA")
  {
    std<-sqrt(diag(beta.i.var))
    
    if(length(X[1,]) > length(x.ols[1,])) optima<-matrix(data=0, nrow=ncol(X), ncol=2, dimnames = list(c(c("Ya",as.character(levels(fixed.fact))),if(is.null(dim(fixed.cov))) deparse(substitute(fixed.cov)) else colnames(fixed.cov), if(is.null(dim(random.cov))) deparse(substitute(random.cov)) else colnames(random.cov)), c("Estimates", "Std.error")))
    else
      
      optima<-matrix(data=0, nrow=ncol(X), ncol=2, dimnames = list(c(as.character(levels(fixed.fact)), if(is.null(dim(fixed.cov))) deparse(substitute(fixed.cov)) else colnames(fixed.cov),if(is.null(dim(random.cov))) deparse(substitute(random.cov)) else colnames(random.cov)), c("Estimates", "Std.error")));
    
    corrected_beta_values<-matrix(data= c(corrected_betas), nrow=(nrow(beta1)), ncol=1, dimnames = list(c(as.character(levels(fixed.fact)), if(is.null(dim(fixed.cov))) deparse(substitute(fixed.cov)) else colnames(fixed.cov),if(is.null(dim(random.cov))) deparse(substitute(random.cov)) else colnames(random.cov)), c("Bias-corr. regression parameters")))
    
    
    optima[,1] = gls.beta1
    optima[,2] = std
    print(optima);message("");message("")
    print(corrected_beta_values);message("")
  }
  if(model.type=="rReg")
  {
    if(ultrametric==TRUE || alpha.est == Inf)
      opreg<-matrix(data=0, nrow=(nrow(gls.beta1)), ncol=2, dimnames=list(c("K", if(is.null(dim(random.cov))) deparse(substitute(random.cov)) else colnames(random.cov)), c("Estimate", "Std. Error")))
    else
    {
      if(alpha.est != Inf)
        opreg<-matrix(data=0, nrow=(nrow(gls.beta1)), ncol=2, dimnames=list(c("Xa", "Bo","Ya" ,if(is.null(dim(random.cov))) deparse(substitute(random.cov)) else colnames(random.cov)), c("Estimate", "Std. Error")))
      else opreg<-matrix(data=0, nrow=(nrow(gls.beta1)), ncol=2, dimnames=list(c("K", if(is.null(dim(random.cov))) deparse(substitute(random.cov)) else colnames(random.cov)), c("Estimate", "Std. Error")))}
    
    opreg[,1] =round(gls.beta1, 5)
    opreg[,2]= round(sqrt(diag(beta.i.var)),5)
    
    
    corrected_beta_values<-matrix(data= c(corrected_betas), nrow=(nrow(beta1)), ncol=1, dimnames=list(c("K", if(is.null(dim(random.cov))) deparse(substitute(random.cov)) else colnames(random.cov)), c("Bias-corr. regression parameters")))
    
    if(model.type=="rReg")
    {
      
      evreg<-matrix(data=0, nrow=(nrow(glsyx.beta1)), ncol=2, dimnames=list(c("Intercept", if(is.null(dim(random.cov))) deparse(substitute(random.cov)) else colnames(random.cov)), c("Estimate", "Std. Error")))
      
      
      evreg[,1] =round(glsyx.beta1, 5)
      evreg[,2]= round(sqrt(diag(ev.beta.i.var)),5)
      
      message("Evolutionary regression"); message("")
      print(evreg);
      message("");
    }
    message("Optimal regression"); message("")
    print(opreg);message("");message("");
    
    print(corrected_beta_values); message("");
    
    if(model.type=="rReg" && ultrametric==TRUE && alpha.est != Inf)
    {
      message("")
      message("Decomposition of K assuming Ya = Xa to get the optimal regression intercept Bo")
      message("")
      
      bo<-opreg[1,1] + (c.factor-1)*(sum(gls.beta1[-1]*theta.X))
      print(bo)
      message("")
      message("(Use this as the intercept when plotting the regression line)")
      
      message("")
    }
  }
  
  
  if(model.type=="mfReg")
  {
    if(ultrametric==TRUE || alpha.est == Inf)
      opreg<-matrix(data=0, nrow=(nrow(gls.beta1)), ncol=2, dimnames=list(c("K",if(is.null(dim(fixed.cov))) deparse(substitute(fixed.cov)) else colnames(fixed.cov),if(is.null(dim(random.cov))) deparse(substitute(random.cov)) else colnames(random.cov)), c("Estimate", "Std. Error")))
    else
    {
      if(alpha.est != Inf)
        opreg<-matrix(data=0, nrow=(nrow(gls.beta1)), ncol=2, dimnames=list(c("Xa", "Bo","Ya" ,if(is.null(dim(fixed.cov))) deparse(substitute(fixed.cov)) else colnames(fixed.cov),if(is.null(dim(random.cov))) deparse(substitute(random.cov)) else colnames(random.cov)), c("Estimate", "Std. Error")))
      else opreg<-matrix(data=0, nrow=(nrow(gls.beta1)), ncol=2, dimnames=list(c("K", if(is.null(dim(fixed.cov))) deparse(substitute(fixed.cov)) else colnames(fixed.cov),if(is.null(dim(random.cov))) deparse(substitute(random.cov)) else colnames(random.cov)), c("Estimate", "Std. Error")))}
    
    opreg[,1] =round(gls.beta1, 5)
    opreg[,2]= round(sqrt(diag(beta.i.var)),5)
    
    if(model.type=="mfReg")
    {
      
      evreg<-matrix(data=0, nrow=(nrow(glsyx.beta1)), ncol=2, dimnames=list(c("Intercept",if(is.null(dim(fixed.cov))) deparse(substitute(fixed.cov)) else colnames(fixed.cov), if(is.null(dim(random.cov))) deparse(substitute(random.cov)) else colnames(random.cov)), c("Estimate", "Std. Error")))
      
      corrected_beta_values<-matrix(data= c(corrected_betas), nrow=(nrow(beta1)), ncol=1, dimnames=list(c("Intercept",if(is.null(dim(fixed.cov))) deparse(substitute(fixed.cov)) else colnames(fixed.cov), if(is.null(dim(random.cov))) deparse(substitute(random.cov)) else colnames(random.cov)), c("Bias-corr. regression parameters")))
      
      evreg[,1] =round(glsyx.beta1, 5)
      evreg[,2]= round(sqrt(diag(ev.beta.i.var)),5)
      
      message("Evolutionary regression"); message("")
      print(evreg);
      message("");
    }
    message("Optimal regression"); message("")
    print(opreg);message("");message("");
    
    print(corrected_beta_values); message("");
    
    if(model.type=="mfReg" && ultrametric==TRUE && alpha.est != Inf)
    {
      message("")
      message("Decomposition of K assuming Ya = Xa to get the optimal regression intercept Bo")
      message("")
      
      bo<-opreg[1,1] + (c.factor-1)*(sum(gls.beta1[-(1:(1+n.fixed.pred))]*theta.X))
      print(bo)
      message("")
      message("(Use this as the intercept when plotting the regression line)")
      
      message("")
    }
  }
  
  message("--------------------------------------------------");
  message("MODEL FIT");message("");
  print(modfit); message("");
  message("==================================================");
  
  print("debug: model.fit.dev")
  
  
} # END OF MODEL FITTING FUNCTION


