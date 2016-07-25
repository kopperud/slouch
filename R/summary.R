## Output

#' Title
#'
#' @param model
#'
#' @return
#' @export
summary.slouch <- function(model = slouch, ...){

  model.type <- model$model.type
  alpha.est <- model$alpha.est
  vy.est <- model$vy.est
  pred.mean <- model$pred.mean
  g.mean <- model$g.mean
  sst <- model$sst
  sse <- model$sse
  r.squared <- model$r.squared
  ml <- model$ml
  n.pred <- model$n.pred
  ultrametric <- model$ultrametric
  intercept <- model$intercept
  beta.est <- model$beta.est
  N <- model$N
  theta.X <- model$theta.X
  s.X <- model$s.X
  beta.var.est <- model$beta.var.est
  X <- model$X
  x.ols <- model$x.ols
  random.cov <- model$random.cov
  corrected_betas <- model$corrected_betas
  fixed.fact <- model$fixed.fact
  random.cov.names <- model$random.cov.names
  glsyx.beta1 <- model$glsyx.beta1
  glsyx.beta1.var <- model$glsyx.beta1.var
  fixed.cov <- model$fixed.cov
  me.fixed.cov <- model$me.fixed.cov
  fixed.cov.names <- model$fixed.cov.names


  gls.beta0 <- gls.beta1 <- beta1 <- beta.est
  beta.i.var <- beta.var.est

  ev.beta.i.var <- glsyx.beta1.var

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
    print(matrix(data=rbind(theta.X, s.X), nrow=2, ncol=n.pred, dimnames=list(c("Predictor theta", "Predictor variance"), random.cov.names)));
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

    optima<-matrix(data=0, nrow=(nrow(gls.beta0)), ncol=2, dimnames=list(c("Bo", fixed.cov.names), c("Estimate", "Std. Error")))
    optima[,1] = gls.beta0;
    optima[,2] = std;

    corrected_beta_values<-matrix(data= c(corrected_betas), nrow=(nrow(beta1)), ncol=1, dimnames=list(c("Bo", fixed.cov.names), c("Bias-corr. regression parameters")))


    print(optima);message("");message("");
    print(corrected_beta_values);message("");

  }

  if(model.type=="ffANCOVA")
  {
    std<-sqrt(diag(beta.i.var))
    optima<-matrix(data=0, nrow=ncol(X), ncol=2, dimnames = list(c(as.character(levels(fixed.fact)), fixed.cov.names), c("Estimates", "Std.error")));
    optima[,1] = gls.beta0;
    optima[,2] = std;

    corrected_beta_values<-matrix(data= c(corrected_betas), nrow=(nrow(beta1)), ncol=1, dimnames = list(c(as.character(levels(fixed.fact)), fixed.cov.names), c("Bias-corr. regression parameters")));

    print(optima);message("");message("");
    print(corrected_beta_values);message("");

  }



  if(model.type  == "mmANCOVA")
  {
    std<-sqrt(diag(beta.i.var))
    if(length(X[1,]) > length(x.ols[1,])) optima<-matrix(data=0, nrow=ncol(X), ncol=2, dimnames = list(c(c("Ya",as.character(levels(fixed.fact))), random.cov.names), c("Estimates", "Std.error")))
    else
      optima<-matrix(data=0, nrow=ncol(X), ncol=2, dimnames = list(c(as.character(levels(fixed.fact)), random.cov.names), c("Estimates", "Std.error")))

    corrected_beta_values<-matrix(data= c(corrected_betas), nrow=(nrow(beta1)), ncol=1, dimnames = list(c(as.character(levels(fixed.fact)), random.cov.names), c("Bias-corr. regression parameters")))


    optima[,1] = gls.beta1;
    optima[,2] = std;
    print(optima);message("");message("");
    print(corrected_beta_values);message("");
  }

  if(model.type  == "mmfANCOVA")
  {
    std<-sqrt(diag(beta.i.var))

    if(length(X[1,]) > length(x.ols[1,])) optima<-matrix(data=0, nrow=ncol(X), ncol=2, dimnames = list(c(c("Ya",as.character(levels(fixed.fact))),fixed.cov.names, random.cov.names), c("Estimates", "Std.error")))
    else

      optima<-matrix(data=0, nrow=ncol(X), ncol=2, dimnames = list(c(as.character(levels(fixed.fact)), fixed.cov.names,random.cov.names), c("Estimates", "Std.error")));

    corrected_beta_values<-matrix(data= c(corrected_betas), nrow=(nrow(beta1)), ncol=1, dimnames = list(c(as.character(levels(fixed.fact)), fixed.cov.names,random.cov.names), c("Bias-corr. regression parameters")))


    optima[,1] = gls.beta1
    optima[,2] = std
    print(optima);message("");message("")
    print(corrected_beta_values);message("")
  }
  if(model.type=="rReg")
  {
    if(ultrametric==TRUE || alpha.est == Inf)
      opreg<-matrix(data=0, nrow=(nrow(gls.beta1)), ncol=2, dimnames=list(c("K", random.cov.names), c("Estimate", "Std. Error")))
    else
    {
      if(alpha.est != Inf)
        opreg<-matrix(data=0, nrow=(nrow(gls.beta1)), ncol=2, dimnames=list(c("Xa", "Bo","Ya" ,random.cov.names), c("Estimate", "Std. Error")))
      else opreg<-matrix(data=0, nrow=(nrow(gls.beta1)), ncol=2, dimnames=list(c("K", random.cov.names), c("Estimate", "Std. Error")))}

    opreg[,1] =round(gls.beta1, 5)
    opreg[,2]= round(sqrt(diag(beta.i.var)),5)

    corrected_beta_values<-matrix(data= c(corrected_betas), nrow=(nrow(beta1)), ncol=1, dimnames=list(c("K", random.cov.names), c("Bias-corr. regression parameters")))

    if(model.type=="rReg")
    {

      evreg<-matrix(data=0, nrow=(nrow(glsyx.beta1)), ncol=2, dimnames=list(c("Intercept", random.cov.names), c("Estimate", "Std. Error")))


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
      opreg<-matrix(data=0, nrow=(nrow(gls.beta1)), ncol=2, dimnames=list(c("K",fixed.cov.names,random.cov.names), c("Estimate", "Std. Error")))
    else
    {
      if(alpha.est != Inf)
        opreg<-matrix(data=0, nrow=(nrow(gls.beta1)), ncol=2, dimnames=list(c("Xa", "Bo","Ya" ,fixed.cov.names,random.cov.names), c("Estimate", "Std. Error")))
      else opreg<-matrix(data=0, nrow=(nrow(gls.beta1)), ncol=2, dimnames=list(c("K", fixed.cov.names,random.cov.names), c("Estimate", "Std. Error")))}

    opreg[,1] =round(gls.beta1, 5)
    opreg[,2]= round(sqrt(diag(beta.i.var)),5)

    if(model.type=="mfReg")
    {
      evreg<-matrix(data=0, nrow=(nrow(glsyx.beta1)), ncol=2, dimnames=list(c("Intercept",fixed.cov.names, random.cov.names), c("Estimate", "Std. Error")))

      corrected_beta_values<-matrix(data= c(corrected_betas), nrow=(nrow(beta1)), ncol=1, dimnames=list(c("Intercept",fixed.cov.names, random.cov.names), c("Bias-corr. regression parameters")))

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


}
