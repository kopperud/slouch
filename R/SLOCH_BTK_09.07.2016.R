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

model.fit<-function(topology, times, half_life_values, vy_values, response, me.response=NULL, fixed.fact=NULL,fixed.cov=NULL, me.fixed.cov=NULL, mecov.fixed.cov=NULL, random.cov=NULL, me.random.cov=NULL, mecov.random.cov=NULL,  intercept="root", ultrametric=TRUE, support=NULL, convergence=NULL, plot.angle=30, parallel.compute = FALSE)
{

  # SET DEFAULTS IF NOT SPECIFIED

  if(is.null(support)) support=2;
  if(is.null(convergence)) convergence=0.000001;
  if(is.null(me.response)) me.response<-diag(rep(0, times=length(response[!is.na(response)])))  else me.response<-diag(me.response[!is.na(me.response)]);


  # DETERMINE MODEL STRUCTURE FROM INPUT AND WRITE A SUMMARY TO THE R CONSOLE

  if(is.null(fixed.fact) && is.null(fixed.cov) && is.null(random.cov)) model.type <- "IntcptReg";
  if(!is.null(fixed.fact) && is.null(fixed.cov) && is.null(random.cov)) model.type <- "ffANOVA";
  if(!is.null(fixed.fact) && !is.null(fixed.cov) && is.null(random.cov)) model.type <-"ffANCOVA";
  if(!is.null(fixed.fact) && is.null(fixed.cov) && !is.null(random.cov)) model.type <- "mmANCOVA";
  if(!is.null(fixed.fact) && !is.null(fixed.cov) && !is.null(random.cov)) model.type <- "mmfANCOVA";
  if(is.null(fixed.fact) && is.null(fixed.cov) && !is.null(random.cov)) model.type <- "rReg";
  if(is.null(fixed.fact) && !is.null(fixed.cov) && is.null(random.cov)) model.type <- "fReg";
  if(is.null(fixed.fact) && !is.null(fixed.cov) && !is.null(random.cov)) model.type <- "mfReg";

  print(model.type)

  # Write type of model to screen

  message("")
  message("MODEL SUMMARY")
  message("")
  if(model.type=="IntcptReg")
  {
    message("You have specified an OU model for a response variable regressed on a grand mean, i.e. one global optima");
    if(ultrametric==FALSE)
    {
      GS_head<-c("Ya", "Theta_Global");
      n.par<-2;
    }
    else
    {
      GS_head<-("Theta_Global");
      n.par<-1;
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

  Y <- response[!is.na(response)];
  N <- length(Y);
  T <- times[terminal.twigs(topology)];
  tia<-tsia(ancestor, time);
  tja<-tsja(ancestor, time);
  term<-terminal.twigs(topology);
  pt<-parse.tree(topology, times);
  ta<-pt$bt;
  tij<-pt$dm;
  num.prob<-matrix(data=0, nrow=N, ncol=N) #this matrix is included for cases where species split at the root;
  cm2<-matrix(data=0, nrow=N, ncol=N);
  gof<-matrix(data=0, nrow=length(half_life_values), ncol=length(vy_values), dimnames=list(half_life_values, vy_values));
  h.lives<-matrix(data=0, nrow=length(half_life_values), ncol=length(vy_values))
  ln2<-log(2)
  half_life_values<-rev(half_life_values)


  # EVALUATE IF IT IS A FIXED FACTOR PREDICTOR OR INTERCEPT ONLY MODEL THEN SET UP APPROPRIATE DESIGN AND VARIANCE MATRICES AND ESTIMATE PARAMETERS WITHOUT ITERATED GLS

  if(model.type =="IntcptReg" || model.type == "ffANOVA")
  {
    if(model.type=="IntcptReg") regime.specs<-rep(1, times=length(topology)) else regime.specs<-fixed.fact;

    cat(c("   ", "t1/2  ", "Vy    ", "Supp  ", GS_head), sep="   ");
    message(" ");


    abc <- cbind(sort(rep(half_life_values, length(vy_values)), decreasing = TRUE), rep(vy_values, length(half_life_values)))
    grid <- lapply(seq_len(nrow(abc)), function(i) abc[i,])
    sup2 <- sapply(grid, hl.intreg, N=N, me.response = me.response, ta = ta, tij = tij, T = T, topology = topology, times = times, regime.specs = regime.specs, model.type = "IntcptReg", ultrametric = TRUE, Y = Y)
    gof <- matrix(sup2, ncol=length(vy_values), byrow=TRUE)

    # print(cbind(abc, sup2))

    # print(grid)

    # Search GOF matrix for best estimates of alpha and vy #

    x<-rev(half_life_values)
    # x <- half_life_values
    y<-vy_values
    z<-gof;
    ml<-max(z);
    for(i in 1:length(half_life_values))
    {
      for(j in 1:length(vy_values))
      {
        if(gof[i,j]==ml){alpha.est=log(2)/half_life_values[i]; vy.est=vy_values[j]}
      }
    }
    for(i in 1:length(half_life_values))
    {
      for(j in 1:length(vy_values))
      {
        if(gof[i,j]<=ml-support) gof[i, j]=ml-support;
      }
    }

    gof=gof-ml

    # final GLS estimations for corrected optima using best alpha and vy estimates #

    if(alpha.est==Inf) alpha.est<-1000000000000000000000

    if(model.type=="IntcptReg")
    {
      if(alpha.est==Inf || alpha.est>=1000000000000000000000 ) X<-matrix(data=1, nrow=N, ncol=1)
      else
        if(ultrametric==TRUE) X<-matrix(data=1, nrow=N, ncol=1)
        else
        {
          X<-matrix(data=0, nrow=N, ncol=2);
          X[,1]<-1-exp(-alpha.est*T);
          X[,2]<-exp(-alpha.est*T)
        }
    }
    else
      X<-weight.matrix(alpha.est, topology,times, N, regime.specs, fixed.cov, intercept)

    V<-((vy.est)*(1-exp(-2*alpha.est*ta))*exp(-alpha.est*tij)) + me.response;
    V.inverse<-solve(V);
    beta.i.var<-pseudoinverse(t(X)%*%V.inverse%*%X);
    beta.i<-beta.i.var%*%(t(X)%*%V.inverse%*%Y);
    gls.beta0<-beta.i;

    # code for calculating SSE, SST and r squared

    pred.mean <- X%*%gls.beta0
    g.mean <- (t(rep(1, times=N))%*%solve(V)%*%Y)/sum(solve(V));
    sst <- t(Y-g.mean)%*% solve(V)%*%(Y-g.mean)
    sse <-t (Y-pred.mean)%*%solve(V)%*%(Y-pred.mean)
    r.squared <- (sst-sse)/sst




  } # END OF FIXED PREDICTOR OR INTERCEPT ONLY PARAMETER ESTIMATION


  if(model.type =="ffANCOVA" || model.type == "fReg")
  {

    fixed.pred<-data.frame(fixed.cov);
    n.fixed.pred<-length(fixed.pred[1,]);
    fixed.pred<-matrix(data=fixed.pred[!is.na(fixed.pred)], ncol=n.fixed.pred);
    if(is.null(me.fixed.cov)) me.fixed.pred<-matrix(data=0, nrow=N, ncol=n.fixed.pred) else me.fixed.pred<- matrix(data=me.fixed.cov[!is.na(me.fixed.cov)], ncol=n.fixed.pred);
    if(is.null(mecov.fixed.cov)) me.cov<-matrix(data=0, nrow=N, ncol=n.fixed.pred) else me.cov<-matrix(data=me.cov.fixed.cov[!is.na(me.cov.fixed.cov)], ncol=n.fixed.pred);


    if(model.type=="fReg")
    {
      x.ols<-cbind(1, fixed.pred);

      beta1<-solve(t(x.ols)%*%x.ols)%*%(t(x.ols)%*%Y);

      n.fixed<-1
      cat(c("   ", "t1/2  ", "Vy    ", "Supp  ", "Bo", if(is.null(dim(fixed.cov))) deparse(substitute(fixed.cov)) else colnames(fixed.cov)), sep="   ");
      message("");
    }

    if(model.type=="ffANCOVA")
    {
      regime.specs<-fixed.fact;
      n.fixed<-length(levels(as.factor(regime.specs)))
      regime.specs<-as.factor(regime.specs)

      x.ols<-weight.matrix(10, topology, times, N, regime.specs, fixed.pred, intercept);


      for (i in length(x.ols[1,]):1){
        if(sum(x.ols[,i]) == 0) {x.ols<-x.ols[,-i]; n.fixed<-n.fixed-1}  #removes internal regimes that have only zero entries

      }
      print("kjetil3")
      beta1<-solve(t(x.ols)%*%x.ols)%*%(t(x.ols)%*%Y);
      print("kjetil4")
      cat(c("   ", "t1/2  ", "Vy    ", "Supp  ", GS_head, if(is.null(dim(fixed.cov))) deparse(substitute(fixed.cov)) else colnames(fixed.cov)), sep="   ");
      message("");
    }

    ## Setting up the Vu and Vd matrices ##

    if(model.type=="fReg") Vu<-diag(c(rep(0,N), c(as.numeric(na.exclude(me.fixed.pred))))) else Vu<-diag(c(rep(0,n.fixed*N),as.numeric(na.exclude(me.fixed.pred))))


    true_var<-matrix(data=0, ncol=n.fixed.pred, nrow=N);
    for (i in 1:n.fixed.pred){
      true_var[,i]<-var(na.exclude(fixed.pred[,i]))-as.numeric(na.exclude(me.fixed.pred[,i]))
    }
    true_var<-c(true_var)


    if(model.type=="fReg") Vd<-diag(c(rep(0,N),true_var)) else Vd<-diag(c(rep(0,n.fixed*N), true_var))

    error_condition<-Vu-(Vu%*%pseudoinverse(Vu+Vd)%*%Vu)

    ## Multiplies with betas ##

    xx<-seq(from=1, to=length(Vu[,1]), by=N)

    obs_var_con <-matrix(0, nrow=N, ncol=N)

    ##
    if(length(xx) > length(beta1)) {beta1<-as.matrix(c(0, beta1))}; #n.fixed<-n.fixed+1 ;}
    #if(length(X[1,])< length(beta1)) {beta1<-solve(t(x.ols)%*%x.ols)%*%(t(x.ols)%*%Y);n.fixed<-length(levels(as.factor(regime.specs))); print("The Ya parameter is dropped as its coefficient is too small");}

    ##

    for (e in seq(from=1, to=ncol(x.ols), by=1)){
      for (j in seq(from=1, to=ncol(x.ols), by=1)) {
        tmp<-error_condition[xx[e]:(e*N),xx[j]:(j*N)]*beta1[e]*beta1[j]
        obs_var_con <-obs_var_con + tmp
      }

    }

    for(i in 1:length(half_life_values))
    {
      for(k in 1:length(vy_values))
      {
        vy <- vy_values[k];

        if(half_life_values[i]==0)
        {
          obs_var_con<-obs_var_con;
        }

        else

          # updates the V matrix by multiplying with new betas #

        {
          obs_var_con <-matrix(0, nrow=N, ncol=N)

          for (e in seq(from=1, to=ncol(x.ols), by=1)){
            for (j in seq(from=1, to=ncol(x.ols), by=1)) {
              tmp<-error_condition[xx[e]:(e*N),xx[j]:(j*N)]*beta1[e]*beta1[j]
              obs_var_con <-obs_var_con + tmp
            }

          }
        }




        if(half_life_values[i]==0)

        {
          a<-1000000000000000000000

          #V<-diag(rep(vy, times=N)) + na.exclude(me.response) + diag(as.numeric(me.fixed.pred%*%(beta1[(n.fixed+1):length(beta1),]*beta1[(n.fixed+1):length(beta1),])))-diag(as.numeric(me.cov%*%(2*beta1[(n.fixed+1):length(beta1),])));


          V<-diag(rep(vy, times=N)) + na.exclude(me.response) + obs_var_con-diag(as.vector(na.exclude(me.cov%*%(2*beta1[(n.fixed+1):length(beta1),]))));



        }

        else

        {

          a <- ln2/half_life_values[i];

          #V<-((vy)*(1-exp(-2*a*ta))*exp(-a*tij))+me.response + diag(as.numeric(me.fixed.pred%*%(beta1[(n.fixed+1):length(beta1),]*beta1[(n.fixed+1):length(beta1),])))-diag(as.numeric(me.cov%*%(2*beta1[(n.fixed+1):length(beta1),])));

          V<-((vy)*(1-exp(-2*a*ta))*exp(-a*tij))+ na.exclude(me.response) + obs_var_con-diag(as.vector(na.exclude(me.cov%*%(2*beta1[(n.fixed+1):length(beta1),]))));

        }

        if(model.type=="fReg") X<-cbind(1, fixed.pred) else  X<-weight.matrix(a, topology, times, N, regime.specs, fixed.pred, intercept);


        ##### iterated GLS



        con.count<-0;  # Counter for loop break if Beta's dont converge #
        repeat
        {

          if(half_life_values[i]==0)
          {
            a<-1000000000000000000000
            #V<-diag(rep(vy, times=N)) + me.response + diag(as.numeric(me.fixed.pred%*%(beta1[(n.fixed+1):length(beta1),]*beta1[(n.fixed+1):length(beta1),])))-diag(as.numeric(me.cov%*%(2*beta1[(n.fixed+1):length(beta1),])));

            V<-((vy)*(1-exp(-2*a*ta))*exp(-a*tij))+ na.exclude(me.response) + obs_var_con-diag(as.vector(na.exclude(me.cov%*%(2*beta1[(n.fixed+1):length(beta1),]))));

          }
          else
          {

            a <- ln2/half_life_values[i];


            obs_var_con <-matrix(0, nrow=N, ncol=N)

            for (e in seq(from=1, to=ncol(x.ols), by=1)){
              for (j in seq(from=1, to=ncol(x.ols), by=1)) {
                tmp<-error_condition[xx[e]:(e*N),xx[j]:(j*N)]*beta1[e]*beta1[j]
                obs_var_con <-obs_var_con + tmp
              }

            }


            #V<-((vy)*(1-exp(-2*a*ta))*exp(-a*tij))+me.response + diag(as.numeric(me.fixed.pred%*%(beta1[(n.fixed+1):length(beta1),]*beta1[(n.fixed+1):length(beta1),])))-diag(as.numeric(me.cov%*%(2*beta1[(n.fixed+1):length(beta1),])));




            V<-((vy)*(1-exp(-2*a*ta))*exp(-a*tij))+ na.exclude(me.response) + obs_var_con-diag(as.vector(na.exclude(me.cov%*%(2*beta1[(n.fixed+1):length(beta1),]))));


          } # END OF If ELSE CONDITION FOR HALF-LIFE 0 OR NOT

          if(model.type=="fReg") X<-cbind(1, fixed.pred) else  X<-weight.matrix(a, topology, times, N, regime.specs, fixed.pred, intercept);

          # INTERMEDIATE ESTIMATION OF OPTIMAL REGRESSION #

          V.inverse<-solve(V)

          beta.i<-pseudoinverse(t(X)%*%V.inverse%*%X)%*%(t(X)%*%V.inverse%*%Y)

          if(length(beta.i) > length(beta1)) {beta1<-as.matrix(c(0, beta1))}

          test<-matrix(nrow=(length(beta.i)))
          for(f in 1:(length(beta.i)))
          {
            if(abs(as.numeric(beta.i[f]-beta1[f]))<=convergence) test[f]=0 else test[f]=1

          }


          if(sum(test)==0) break
          con.count=con.count+1
          if(con.count >= 50)
          {
            message("Warning, estimates did not converge after 50 iterations, last estimates printed out")
            break
          }
          beta1<-beta.i


        }


        eY<-X%*%beta1
        resid<-Y-eY

        det.V<-det(V)
        if(det.V==0){
          print(paste("Warning: Determinant of V = 0"))
          #Minimum value of diagonal scaling factor
          inv.min.diag.V<-1/min(diag(V))
          V<-V*inv.min.diag.V
          #Rescale and log determinant
          log.det.V<-log(det(V))+log(min(diag(V)))*N
        }
        else {log.det.V<-log(det.V)}

        gof[i, k] <- -N/2*log(2*pi)-0.5*log.det.V-0.5*(t(resid) %*% V.inverse%*%resid);
        print(as.numeric(round(cbind(if(a!=0)log(2)/a else 0.00, vy, gof[i,k], t(beta1)), 4)))



        ### END OF ITERATED GLS

      } # end of half-life loop
    } # end of vy loop

    # Search GOF matrix for best estimates of alpha and vy #

    x<-rev(half_life_values)
    y<-vy_values
    z<-gof;
    ml<-max(z);
    for(i in 1:length(half_life_values))
    {
      for(j in 1:length(vy_values))
      {
        if(gof[i,j]==ml){alpha.est=log(2)/half_life_values[i]; vy.est=vy_values[j]}
      }
    }
    for(i in 1:length(half_life_values))
    {
      for(j in 1:length(vy_values))
      {
        if(gof[i,j]<=ml-support) gof[i, j]=ml-support;
      }
    }

    gof=gof-ml



    # final GLS estimations for corrected optima using best alpha and vy estimates #



    #x.ols<-cbind(1, fixed.pred);
    #beta1<-solve(t(x.ols)%*%x.ols)%*%(t(x.ols)%*%Y);




    con.count<-0;  # Counter for loop break if Beta's dont converge #
    repeat
    {
      if(alpha.est==Inf)
      {
        a<-1000000000000000000000

        #V<-diag(rep(vy, times=N)) + me.response + diag(as.numeric(me.fixed.pred%*%(beta1[(n.fixed+1):length(beta1),]*beta1[(n.fixed+1):length(beta1),])))-diag(as.numeric(me.cov%*%(2*beta1[(n.fixed+1):length(beta1),])));

        V<-diag(rep(vy.est, times=N)) + na.exclude(me.response) + obs_var_con-diag(as.vector(na.exclude(me.cov%*%(2*beta1[(n.fixed+1):length(beta1),]))));

      }
      else
      {

        #V<-((vy)*(1-exp(-2*alpha.est*ta))*exp(-alpha.est*tij))+me.response + diag(as.numeric(me.fixed.pred%*%(beta1[(n.fixed+1):length(beta1),]*beta1[(n.fixed+1):length(beta1),])))-diag(as.numeric(me.cov%*%(2*beta1[(n.fixed+1):length(beta1),])));



        obs_var_con <-matrix(0, nrow=N, ncol=N)

        for (e in seq(from=1, to=ncol(x.ols), by=1)){
          for (j in seq(from=1, to=ncol(x.ols), by=1)) {
            tmp<-error_condition[xx[e]:(e*N),xx[j]:(j*N)]*beta1[e]*beta1[j]
            obs_var_con <-obs_var_con + tmp
          }

        }

        V<-((vy.est)*(1-exp(-2*alpha.est*ta))*exp(-alpha.est*tij))+  na.exclude(me.response) + obs_var_con-diag(as.vector(na.exclude(me.cov%*%(2*beta1[(n.fixed+1):length(beta1),]))));


      } # END OF If ELSE CONDITION FOR HALF-LIFE 0 OR NOT

      if(model.type=="fReg") X<-cbind(1, fixed.pred) else  X<-weight.matrix(a, topology, times, N, regime.specs, fixed.pred, intercept);

      # INTERMEDIATE ESTIMATION OF OPTIMAL REGRESSION #

      V.inverse<-solve(V)

      beta.i<-pseudoinverse(t(X)%*%V.inverse%*%X)%*%(t(X)%*%V.inverse%*%Y)
      test<-matrix(nrow=(length(beta.i)))
      for(f in 1:(length(beta.i)))
      {
        if(abs(as.numeric(beta.i[f]-beta1[f]))<=convergence) test[f]=0 else test[f]=1
      }
      if(sum(test)==0) break
      con.count=con.count+1
      if(con.count >= 50)
      {
        message("Warning, estimates did not converge after 50 iterations, last estimates printed out")
        break
      }

      beta1<-beta.i
    }


    gls.beta0<-beta1;

    beta.i.var<-pseudoinverse(t(X)%*%V.inverse%*%X);

    # code for calculating SSE, SST and r squared

    pred.mean <- X%*%gls.beta0
    g.mean <- (t(rep(1, times=N))%*%solve(V)%*%Y)/sum(solve(V));
    sst <- t(Y-g.mean)%*% solve(V)%*%(Y-g.mean)
    sse <-t (Y-pred.mean)%*%solve(V)%*%(Y-pred.mean)
    r.squared <- (sst-sse)/sst



    ###### Bias correction ######

    if(model.type=="fReg") X<-cbind(1, fixed.pred) else  X<-weight.matrix(alpha.est, topology, times, N, regime.specs, fixed.pred, intercept);

    adj<-matrix(data=0, ncol=ncol(X), nrow=N)  #PREDICTOR THETA
    for(i in 1:length(adj[1,]))
    {
      adj[,i] <- mean(X[,i]);
    }

    correction<-matrix(Vu%*%pseudoinverse(Vd+Vu)%*%(c(X)-c(adj)),  ncol=ncol(X), nrow=nrow(X), byrow=F)

    bias_corr<-pseudoinverse(t(X)%*%V.inverse%*%X)%*%t(X)%*%V.inverse%*%correction


    m<-length(gls.beta0);



    corrected_betas<-solve(diag(1,m,m)-bias_corr)%*%gls.beta0

    ###### End of Bias correction ######


  }  # END OF fReg AND ffANVOCA ESTIMATION ROUTINES must still add iterated GLS for me



  # EVALUATE IF IT IS A FIXED MODEL ANCOVA, MIXED MODEL ANCOVA OR RANDOM PREDICTOR REGRESSION, ESTIMATE PARAMETERS WITH ITERATED GLS TO A) TAKE MEASUREMENT VARIANCE INTO ACCOUNT OR B) RANDOM EFFECTS INTO ACCOUNT IN THE CASE OF THE MIXED MODEL AND REGRESSION

  if(model.type == "mmANCOVA" || model.type=="rReg")  ### more models here
  {
    # SET UP INITIAL MATRICES FOR MULTIPLE REGRESSION AND CALCULATE THETA AND SIGMA FOR RANDOM PREDICTOR / S

    pred<-data.frame(random.cov);
    n.pred<-length(pred[1,]);
    pred<-matrix(data=pred[!is.na(pred)], ncol=n.pred);
    if(is.null(me.random.cov)) me.pred<-matrix(data=0, nrow=N, ncol=n.pred) else me.pred<-matrix(data=me.random.cov[!is.na(me.random.cov)], ncol=n.pred);
    if(is.null(mecov.random.cov)) me.cov<-matrix(data=0, nrow=N, ncol=n.pred) else me.cov<-matrix(data=mecov.random.cov[!is.na(mecov.random.cov)], ncol=n.pred);

    s.X<-matrix(data=0, ncol=n.pred)  # PREDICTOR SIGMA
    for(i in 1:n.pred)
    {
      s.X[,i] <- as.numeric(sigma.X.estimate(pred[,i], me.pred[,i], topology, times)[2]);
    }

    theta.X<-matrix(data=0, ncol=n.pred)  #PREDICTOR THETA
    for(i in 1:n.pred)
    {
      theta.X[,i] <- as.numeric(sigma.X.estimate(pred[,i],me.pred[,i], topology, times)[1]);
    }



    # END OF RANDOM PREDICTOR THETA AND SIGMA ESTIMATES


    ## INITIAL OLS ESTIMATES TO SEED ITERATED GLS

    if(model.type=="rReg")
    {
      x.ols<-cbind(1, pred);
      beta1<-solve(t(x.ols)%*%x.ols)%*%(t(x.ols)%*%Y);
      if(ultrametric == FALSE) beta1<-rbind(0, 0, beta1); # 2 additional parameter seeds for Ya and Xa
      n.fixed<-1

      ## Setting up the Vu and Vd matrices ##


      Vd<-matrix(0,ncol=(N*length(beta1[,1])), nrow=(N*length(beta1[,1])))

      xx<-seq(from=1, to	=length(Vd[,1]), by=N)

      if(ultrametric == TRUE) xx<-xx[-1] else xx<-xx[-(1:3)]

      yy<-seq(from=N, to	=length(Vd[,1]), by=N)

      if(ultrametric == TRUE) yy<-yy[-1] else yy<-yy[-(1:3)]


      for (i in seq(from=1, to=nrow(s.X), by=1)){
        Vd[xx[i]:yy[i],xx[i]:yy[i]]<-pt$bt*s.X[,i]

      }
      if(ultrametric == TRUE) Vu<-diag(c(rep(0,N), c(as.numeric(na.exclude(me.pred))))) else Vu<-diag(c(rep(0,N*3), c(as.numeric(na.exclude(me.pred)))))

      error_condition<-Vu-(Vu%*%pseudoinverse(Vu+Vd)%*%Vu)


      xx<-seq(from=1, to=length(Vu[,1]), by=N)

      obs_var_con <-matrix(0, nrow=N, ncol=N)

      for (e in seq(from=1, to=ncol(x.ols), by=1)){
        for (j in seq(from=1, to=ncol(x.ols), by=1)) {
          tmp<-error_condition[xx[e]:(e*N),xx[j]:(j*N)]*beta1[e]*beta1[j]
          obs_var_con <-obs_var_con + tmp
        }

      }


    }




    if(model.type=="mmANCOVA")
    {
      regime.specs<-fixed.fact;
      n.fixed<-length(levels(as.factor(regime.specs)))
      regime.specs<-as.factor(regime.specs)
      x.ols<-cbind(weight.matrix(1000000000000000000000, topology, times, N, regime.specs, fixed.cov=NULL, intercept), pred);

      for (i in length(x.ols[1,]):1){
        if(sum(x.ols[,i]) == 0) {x.ols<-x.ols[,-i]; n.fixed<-n.fixed-i}  #removes internal regimes that have only zero entries

      }
      beta1<-solve(t(x.ols)%*%x.ols)%*%(t(x.ols)%*%Y);

      ## Setting up the Vu and Vd matrices ##


      Vd<-matrix(0,ncol=(N*length(beta1[,1])), nrow=(N*length(beta1[,1])))

      xx<-seq(from=1, to	=length(Vd[,1]), by=N)

      #if(ultrametric == TRUE) xx<-xx[-1] else

      xx<-xx[-(1: n.fixed)]

      yy<-seq(from=N, to	=length(Vd[,1]), by=N)

      #if(ultrametric == TRUE) yy<-yy[-1] else

      yy<-yy[-(1: n.fixed)]


      for (i in seq(from=1, to=nrow(s.X), by=1)){
        Vd[xx[i]:yy[i],xx[i]:yy[i]]<-pt$bt*s.X[,i]
      }

      #if(ultrametric == TRUE) Vu<-diag(c(rep(0,N), c(as.numeric(na.exclude(me.pred))))) else

      Vu<-diag(c(rep(0,N* n.fixed), c(as.numeric(na.exclude(me.pred)))))

      error_condition<-Vu-(Vu%*%pseudoinverse(Vu+Vd)%*%Vu)


      xx<-seq(from=1, to=length(Vu[,1]), by=N)
      obs_var_con <-matrix(0, nrow=N, ncol=N)

      for (e in seq(from=1, to=ncol(x.ols), by=1)){
        for (j in seq(from=1, to=ncol(x.ols), by=1)) {
          tmp<-error_condition[xx[e]:(e*N),xx[j]:(j*N)]*beta1[e]*beta1[j]
          obs_var_con <-obs_var_con + tmp
        }

      }

    }


    # GRID ESTIMATION ROUTINE AND ITERATED GLS FOR MODELS THAT INCLUDE RANDOM EFFECTS

    if(model.type=="mmANCOVA")
    {

      cat(c("   ", "t1/2  ", "Vy    ", "Supp  ", GS_head, if(is.null(dim(random.cov))) deparse(substitute(random.cov)) else colnames(random.cov)), sep="   ");


      message(" ");


      for(i in 1:length(half_life_values))
      {
        for(k in 1:length(vy_values))
        {



          if(half_life_values[i]==0) a<-1000000000000000000000 else a <- ln2/half_life_values[i];

          if(half_life_values[i]==0)
          {
            obs_var_con<-obs_var_con;
          }

          else

            #multiplying the V matrix with new betas

          {
            obs_var_con <-matrix(0, nrow=N, ncol=N)

            for (e in seq(from=1, to=ncol(x.ols), by=1)){
              for (j in seq(from=1, to=ncol(x.ols), by=1)) {
                tmp<-error_condition[xx[e]:(e*N),xx[j]:(j*N)]*beta1[e]*beta1[j]
                obs_var_con <-obs_var_con + tmp
              }

            }
          }

          vy <- vy_values[k];
          X<-cbind(weight.matrix(a, topology, times, N, regime.specs, fixed.cov=NULL, intercept), (1-(1-exp(-a*T))/(a*T))*pred);
          if(length(X[1,]) > length(beta1)) {beta1<-as.matrix(c(0, beta1)); n.fixed<-n.fixed+1 ;}
          if(length(X[1,])< length(beta1)) {beta1<-solve(t(x.ols)%*%x.ols)%*%(t(x.ols)%*%Y);n.fixed<-length(levels(as.factor(regime.specs))); print("The Ya parameter is dropped as its coefficient is too small");}


          # CODE FOR ESTIMATING BETA USING ITERATED GLS

          con.count<-0;  # Counter for loop break if Beta's dont converge #
          repeat
          {
            if(half_life_values[i]==0)
            {
              X<-cbind(weight.matrix(1000000000000000000000, topology, times, N, regime.specs, fixed.cov=NULL, intercept), pred);
              #V<-diag(rep(vy, times=N))+me.response+diag(as.numeric(me.pred%*%(beta1[(n.fixed+1):length(beta1),]*beta1[(n.fixed+1):length(beta1),])))-diag(as.numeric(me.cov%*%(2*beta1[(n.fixed+1):length(beta1),])));

              V<-diag(rep(vy, times=N))+na.exclude(me.response)+ obs_var_con-diag(as.numeric(me.cov%*%(2*beta1[(n.fixed+1):length(beta1),])));

            }
            else
            {


              X<-cbind(weight.matrix(a, topology, times, N, regime.specs, fixed.cov=NULL, intercept), (1-(1-exp(-a*T))/(a*T))*pred);


              s1<-as.numeric(s.X%*%(beta1[(n.fixed+1):length(beta1),]*beta1[(n.fixed+1):length(beta1),]));


              for(p in 1:N)
              {
                for(q in 1:N)
                {
                  if(ta[q,p]==0)num.prob[q,p]=1 else num.prob[q,p]=(1-exp(-a*ta[q,p]))/(a*ta[q,p]);
                }
              }
              cm1<-(s1/(2*a)+vy)*(1-exp(-2*a*ta))*exp(-a*tij);
              for(p in 1:N)
              {
                for(q in 1:N)
                {
                  cm2[p,q]<-(((1-exp(-a*T[p]))/(a*T[p]))*((1-exp(-a*T[q]))/(a*T[q]))-(exp(-a*tia[p, q])*(1-exp(-a*T[p]))/ (a*T[q])+ exp(-a*tja[p, q])*(1-exp(-a*T[p]))/(a*T[p]))*(num.prob[p,q]));
                }
              }

              mv<-diag(rowSums(matrix(data=as.numeric(me.pred)*t(kronecker(beta1[(n.fixed+1):length(beta1), ], (1-(1-exp(-a*T))/(a*T)))^2), ncol=n.pred)));
              mcov<-diag(rowSums(matrix(data=as.numeric(me.cov)*t(kronecker(2*beta1[(n.fixed+1):length(beta1),], (1-(1-exp(-a*T))/(a*T)))), ncol=n.pred)));


              #V<-cm1+(s1*ta*cm2)+me.response+mv-mcov

              V<-cm1+(s1*ta*cm2)+na.exclude(me.response)+ obs_var_con-mcov


            } # END OF If ELSE CONDITION FOR HALF-LIFE 0 OR NOT

            # INTERMEDIATE ESTIMATION OF OPTIMAL REGRESSION #

            V.inverse<-solve(V)
            if(half_life_values[i]==0)
            {

              for (z in length(X[1,]):1)
              {
                if(sum(X[,z]) == 0) {X<-X[,-z]; n.fixed-z}  #removes internal regimes that have only zero entries
              }


              beta.i<-pseudoinverse(t(X)%*%V.inverse%*%X)%*%(t(X)%*%V.inverse%*%Y)
              test<-matrix(nrow=(length(beta.i)))
              for(f in 1:(length(beta.i)))
              {
                if(abs(as.numeric(beta.i[f]-beta1[f]))<=convergence) test[f]=0 else test[f]=1
              }
              if(sum(test)==0) break
              con.count=con.count+1
              if(con.count >= 50)
              {
                message("Warning, estimates did not converge after 50 iterations, last estimates printed out")
                break
              }

              beta1<-beta.i

            }
            else
            {

              beta.i<-pseudoinverse(t(X)%*%V.inverse%*%X)%*%(t(X)%*%V.inverse%*%Y)
              test<-matrix(nrow=(length(beta.i)))
              for(f in 1:(length(beta.i)))
              {
                if(abs(as.numeric(beta.i[f]-beta1[f]))<=convergence) test[f]=0 else test[f]=1
              }
              if(sum(test)==0) break
              con.count=con.count+1
              if(con.count >= 50)
              {
                message("Warning, estimates did not converge after 50 iterations, last estimates printed out")
                break
              }

              beta1<-beta.i

            }
          }

          ### END OF ITERATED GLS ESTIMATION FOR BETA #
          if(half_life_values[i]==0)
          {
            print(half_life_values[i]);
            X<-cbind(weight.matrix(1000000000000000000000, topology, times, N, regime.specs, fixed.cov=NULL, intercept), pred)
            n.fixed<-length(levels(as.factor(regime.specs)))


            if(length(X[1,]) > length(beta1)) {beta1<-as.matrix(c(0, beta1));}
            if(length(X[1,])< length(beta1)) {beta1<-solve(t(x.ols)%*%x.ols)%*%(t(x.ols)%*%Y);n.fixed<-length(levels(as.factor(regime.specs))); print("The Ya parameter is dropped as its coefficient is too small");}

            #V<-diag(rep(vy, times=N))+me.response+diag(as.numeric(me.pred%*%(beta1[(n.fixed+1):length(beta1),]*beta1[(n.fixed+1):length(beta1),])))-diag(as.numeric(me.cov%*%(2*beta1[(n.fixed+1):length(beta1),])))



            V<-diag(rep(vy, times=N))+na.exclude(me.response)+ obs_var_con-diag(as.numeric(me.cov%*%(2*beta1[(n.fixed+1):length(beta1),])))

            # k<-diag(rep(vy, times=N))+me.response+ obs_var_con-diag(as.numeric(me.cov%*%(2*beta1[(n.fixed+1):length(beta1),])))



            V.inverse<-solve(V)
            eY<-X%*%beta1
            resid<-Y-eY;

            det.V<-det(V)
            if(det.V==0){
              print(paste("Warning: Determinant of V = 0"))
              #Minimum value of diagonal scaling factor
              inv.min.diag.V<-1/min(diag(V))
              V<-V*inv.min.diag.V
              #Rescale and log determinant
              log.det.V<-log(det(V))+log(min(diag(V)))*N
            }
            else {log.det.V<-log(det.V)}

            gof[i, k] <- -N/2*log(2*pi)-0.5*log(det(V))-0.5*(t(resid) %*% V.inverse%*%resid);
          }
          else
          {
            s1<-as.numeric(s.X%*%(beta1[(n.fixed+1):length(beta1),]*beta1[(n.fixed+1):length(beta1),]))
            for(p in 1:N)
            {
              for(q in 1:N)
              {
                if(ta[q,p]==0)num.prob[q,p]=1 else num.prob[q,p]=(1-exp(-a*ta[q,p]))/(a*ta[q,p]);
              }
            }
            cm1<-(s1/(2*a)+vy)*(1-exp(-2*a*ta))*exp(-a*tij);
            for(p in 1:N)
            {
              for(q in 1:N)
              {
                cm2[p,q]<-(((1-exp(-a*T[p]))/(a*T[p]))*((1-exp(-a*T[q]))/(a*T[q]))-(exp(-a*tia[p, q])*(1-exp(-a*T[p]))/(a*T[q])+ exp(-a*tja[p, q])*(1-exp(-a*T[p]))/(a*T[p]))*(num.prob[p,q]));
              }
            }
            X<-cbind(weight.matrix(a, topology, times, N, regime.specs, fixed.cov, intercept), (1-(1-exp(-a*T))/(a*T))*pred);
            mv<-diag(rowSums(matrix(data=as.numeric(me.pred)*t(kronecker(beta1[(n.fixed+1):length(beta1), ], (1-(1-exp(-a*T))/(a*T)))^2), ncol=n.pred)));
            mcov<-diag(rowSums(matrix(data=as.numeric(me.cov)*t(kronecker(2*beta1[(n.fixed+1):length(beta1), ], (1-(1-exp(-a*T))/(a*T)))), ncol=n.pred)));
            #V<-cm1+(s1*ta*cm2)+me.response+mv-mcov;
            V<-cm1+(s1*ta*cm2)+na.exclude(me.response)+ obs_var_con-mcov


            V.inverse<-solve(V)

            eY<-X%*%beta1

            resid<-Y-eY;

            det.V<-det(V)
            if(det.V==0){
              print(paste("Warning: Determinant of V = 0"))
              #Minimum value of diagonal scaling factor
              inv.min.diag.V<-1/min(diag(V))
              V<-V*inv.min.diag.V
              #Rescale and log determinant
              log.det.V<-log(det(V))+log(min(diag(V)))*N
            }
            else {log.det.V<-log(det.V)}

            gof[i, k] <- -N/2*log(2*pi)-0.5*log.det.V-0.5*(t(resid) %*% V.inverse%*%resid);
          }  # END OF CONDITION FOR HALF-LIFE = 0 #
          print(as.numeric(round(cbind(if(a!=0)log(2)/a else 0.00, vy, gof[i,k], t(beta1)), 4)))

        }
      }



      # END OF GRID SETUP,START OF GRID SEARCH FOR BEST ALPHA AND VY ESTIMATES #

      x<-rev(half_life_values)
      y<-vy_values
      z<-gof;
      ml<-max(z);
      for(i in 1:length(half_life_values))
      {
        for(j in 1:length(vy_values))
        {
          if(gof[i,j]==ml){alpha.est=log(2)/half_life_values[i]; vy.est=vy_values[j]}
        }
      }
      for(i in 1:length(half_life_values))
      {
        for(j in 1:length(vy_values))
        {
          if(gof[i,j]<=ml-support)gof[i, j]=ml-support;
        }
      }
      gof=gof-ml

      n.fixed<-length(levels(as.factor(regime.specs)))   ### reset before final regression



      # FINAL OPTIMAL REGRESSION USING BEST ALPHA AND VY ESTIMATES #

      if(alpha.est==Inf || alpha.est >=1000000000000000000000)
      {
        x.ols<-cbind(weight.matrix(1000000000000000000000, topology, times, N, regime.specs, fixed.cov, intercept), pred)
        gls.beta1<-solve(t(x.ols)%*%x.ols)%*%(t(x.ols)%*%Y)
        con.count<-0;
        repeat
        {
          s1<-as.numeric(s.X%*%(gls.beta1[(n.fixed+1):length(gls.beta1),]*gls.beta1[(n.fixed+1):length(gls.beta1),]))
          X<-cbind(weight.matrix(1000000000000000000000, topology, times, N, regime.specs, fixed.cov, intercept), pred)
          #V<-diag(rep(vy, times=N))+me.response+diag(as.numeric(me.pred%*%(gls.beta1[(n.fixed+1):length(gls.beta1),]*gls.beta1[(n.fixed+1):length(gls.beta1),])))-diag(as.numeric(me.cov%*%(2*gls.beta1[(n.fixed+1):length(gls.beta1),])))

          V<-diag(rep(vy.est, times=N))+na.exclude(me.response)+ obs_var_con-diag(as.numeric(me.cov%*%(2*gls.beta1[(n.fixed+1):length(gls.beta1),])))


          V.inverse<-solve(V)
          beta.i.var<-ev.beta.i.var<-pseudoinverse(t(X)%*%V.inverse%*%X)
          beta.i<-beta.i.var%*%(t(X)%*%V.inverse%*%Y)
          test<-matrix(nrow=(length(beta.i)))
          for(f in 1:(length(beta.i)))
          {
            if(abs(as.numeric(beta.i[f]-gls.beta1[f]))<=convergence) test[f]=0 else test[f]=1
          }
          if(sum(test)==0) break
          con.count=con.count+1
          if(con.count >= 50)
          {
            message("Warning, estimates did not converge after 50 iterations, last estimates printed out")
            break
          }
          gls.beta1<-beta.i
        }
        gls.beta1<-beta.i
        X<-cbind(weight.matrix(1000000000000000000000, topology, times, N, regime.specs, fixed.cov, intercept), pred)
        #V<-diag(rep(vy, times=N))+me.response+diag(as.numeric(me.pred%*%(gls.beta1[(n.fixed+1):length(beta1),]*gls.beta1[(n.fixed+1):length(gls.beta1),])))-diag(as.numeric(me.cov%*%(2*gls.beta1[(n.fixed+1):length(gls.beta1),])))

        V<-diag(rep(vy.est, times=N))+na.exclude(me.response)+ obs_var_con-diag(as.numeric(me.cov%*%(2*gls.beta1[(n.fixed+1):length(gls.beta1),])))

        pred.mean<-X%*%gls.beta1
        g.mean<-(t(rep(1, times=N))%*%solve(V)%*%Y)/sum(solve(V));
        sst<-t(Y-g.mean)%*% solve(V)%*%(Y-g.mean)

        sse<-t(Y-pred.mean)%*%solve(V)%*%(Y-pred.mean)

        r.squared<-(sst-sse)/sst


      }
      else
      {
        x.ols<-cbind(weight.matrix(1000000000000000000000, topology, times, N, regime.specs, fixed.cov, intercept), pred)
        gls.beta1<-solve(t(x.ols)%*%x.ols)%*%(t(x.ols)%*%Y)
        con.count<-0;

        X<-cbind(weight.matrix(alpha.est, topology, times, N, regime.specs, fixed.cov, intercept), (1-(1-exp(-alpha.est*T))/(alpha.est*T))*pred);
        if(length(X[1,]) > length(gls.beta1)) {gls.beta1<-as.matrix(c(0, gls.beta1)); n.fixed<-n.fixed+1}
        if(length(X[1,])< length(gls.beta1)) {gls.beta1<-solve(t(x.ols)%*%x.ols)%*%(t(x.ols)%*%Y);n.fixed<-length(levels(as.factor(regime.specs)))}

        repeat
        {

          X<-cbind(weight.matrix(alpha.est, topology, times, N, regime.specs, fixed.cov, intercept), (1-(1-exp(-alpha.est*T))/(alpha.est*T))*pred);
          s1<-as.numeric(s.X%*%(gls.beta1[(n.fixed+1):length(gls.beta1),]*gls.beta1[(n.fixed+1):length(gls.beta1),]))

          obs_var_con <-matrix(0, nrow=N, ncol=N)

          for (e in seq(from=1, to=ncol(x.ols), by=1)){
            for (j in seq(from=1, to=ncol(x.ols), by=1)) {
              tmp<-error_condition[xx[e]:(e*N),xx[j]:(j*N)]* (gls.beta1[e]*(1-(1-exp(-alpha.est*T))/(alpha.est*T)))*(gls.beta1[j]*(1-(1-exp(-alpha.est*T))/(alpha.est*T)))
              obs_var_con <-obs_var_con + tmp
            }

          }


          for(p in 1:N)
          {
            for(q in 1:N)
            {
              if(ta[q,p]==0)num.prob[q,p]=1 else num.prob[q,p]=(1-exp(-alpha.est*ta[q,p]))/(alpha.est*ta[q,p])
            }
          }
          cm1<-(s1/(2*alpha.est)+vy.est)*(1-exp(-2*alpha.est*ta))*exp(-alpha.est*tij)
          for(p in 1:N)
          {
            for(q in 1:N)
            {
              cm2[p,q]<-(((1-exp(-alpha.est*T[p]))/(alpha.est*T[p]))*((1-exp(-alpha.est*T[q]))/(alpha.est*T[q]))-(exp(-alpha.est*tia[p, q])*(1-exp(-alpha.est*T[p]))/(alpha.est*T[q])+ exp(-alpha.est*tja[p, q])*(1-exp(-alpha.est*T[p]))/(alpha.est*T[p]))*(num.prob[p,q]))
            }
          }


          mv<-diag(rowSums(matrix(data=as.numeric(me.pred)*t(kronecker(gls.beta1[(n.fixed+1):length(gls.beta1), ], (1-(1-exp(-alpha.est*T))/(alpha.est*T)))^2), ncol=n.pred)))
          mcov<-diag(rowSums(matrix(data=as.numeric(me.cov)*t(kronecker(gls.beta1[(n.fixed+1):length(gls.beta1),], (1-(1-exp(-alpha.est*T))/(alpha.est*T)))*2), ncol=n.pred)))
          #V<-cm1+(s1*ta*cm2)+me.response+mv-mcov;
          V<-cm1+(s1*ta*cm2)+na.exclude(me.response)+ obs_var_con-mcov

          V.inverse<-solve(V)
          beta.i.var<-pseudoinverse(t(X)%*%V.inverse%*%X)
          beta.i<-beta.i.var%*%(t(X)%*%V.inverse%*%Y)
          test<-matrix(nrow=(length(beta.i)))
          for(f in 1:(length(beta.i)))
          {
            if(abs(as.numeric(beta.i[f]-gls.beta1[f]))<=convergence) test[f]=0 else test[f]=1
          }
          if(sum(test)==0) break
          con.count=con.count+1
          if(con.count >= 50)
          {
            message("Warning, estimates did not converge after 50 iterations, last estimates printed out")
            break
          }


          gls.beta1<-beta.i

          obs_var_con <-matrix(0, nrow=N, ncol=N)

          for (e in seq(from=1, to=ncol(x.ols), by=1)){
            for (j in seq(from=1, to=ncol(x.ols), by=1)) {
              tmp<-error_condition[xx[e]:(e*N),xx[j]:(j*N)]* (gls.beta1[e]*(1-(1-exp(-alpha.est*T))/(alpha.est*T)))*(gls.beta1[j]*(1-(1-exp(-alpha.est*T))/(alpha.est*T)))
              obs_var_con <-obs_var_con + tmp
            }

          }

          X<-cbind(weight.matrix(alpha.est, topology, times, N, regime.specs, fixed.cov=NULL, intercept), (1-(1-exp(-alpha.est*T))/(alpha.est*T))*pred)



          mv<-diag(rowSums(matrix(data=as.numeric(me.pred)*t(kronecker(gls.beta1[(n.fixed+1):length(gls.beta1), ], (1-(1-exp(-alpha.est*T))/(alpha.est*T)))^2), ncol=n.pred)))
          mcov<-diag(rowSums(matrix(data=as.numeric(me.cov)*t(kronecker(gls.beta1[(n.fixed+1):length(gls.beta1), ], (1-(1-exp(-alpha.est*T))/(alpha.est*T)))*2), ncol=n.pred)))
          #V<-cm1+(s1*ta*cm2)+me.response+mv-mcov;
          V<-cm1+(s1*ta*cm2)+na.exclude(me.response)+ obs_var_con-mcov

          pred.mean<-X%*%gls.beta1
          g.mean<-(t(rep(1, times=N))%*%solve(V)%*%Y)/sum(solve(V));
          sst<-t(Y-g.mean)%*% solve(V)%*%(Y-g.mean)
          sse<-t(Y-pred.mean)%*%solve(V)%*%(Y-pred.mean)
          r.squared<-(sst-sse)/sst




        }
      }

      # END OF ITERATED GLS LOOP #


      ###### Start of Bias correction ######


      adj<-matrix(data=0, ncol= ncol(X), nrow=N)



      for(i in 1:n.pred)
      {
        adj[,(n.fixed+i)] <- as.numeric(sigma.X.estimate(pred[,i], me.pred[,i], topology, times)[1]);
      }

      V.inverse<-solve(V)

      correction<-matrix(Vu%*%pseudoinverse(Vd+Vu)%*%(c(X)-c(adj)),  ncol=ncol(X), nrow=nrow(X), byrow=F)



      bias_corr<-pseudoinverse(t(X)%*%V.inverse%*%X)%*%t(X)%*%V.inverse%*%correction


      m<-length(gls.beta1);



      K<-solve(diag(1,m,m)-bias_corr)


      corrected_betas<-solve(diag(1,m,m)-bias_corr)%*% gls.beta1



      ###### End of Bias correction ######




    } # END OF ESTIMATION MIXED MODEL ANCOVA


    if(model.type=="rReg")
    {

      if(ultrametric==TRUE)
      {

        cat(c("   ", "t1/2  ", "Vy    ", "Supp  ", "K     ", if(is.null(dim(random.cov))) deparse(substitute(random.cov)) else colnames(random.cov)), sep="   ");

      }

      else
        cat(c("   ", "t1/2  ", "Vy    ", "Supp  ", "Ya    ", "Xa    " ,"Bo    ",  if(is.null(dim(random.cov))) deparse(substitute(random.cov)) else colnames(random.cov)), sep="   ");

      message(" ");

      grid <- cbind(sort(rep(half_life_values, length(vy_values)), decreasing = TRUE), rep(vy_values, length(half_life_values)))
      # grid2 <- expand.grid(half_life_values, vy_values); grid2 <- grid2[order(-grid2[,1], grid[,2])]

      if (parallel.compute){
        cl.cores.tmp <- detectCores()
        if (!is.na(cl.cores.tmp)) cl.cores <- cl.cores.tmp else cl.cores <- 1
        cl <- makeCluster(getOption("cl.cores",cl.cores))
        # clusterExport(cl, envir = environment(), varlist=c("N", "me.response", "ta", "tij", "T", "topology", "times", "model.type", "ultrametric", "Y", "fixed.cov", "pred", "xx", "beta1", "error_condition", "s.X", "n.pred", "num.prob", "tia", "tja", "cm2", "me.pred", "me.cov", "convergence", "n.fixed"))
        sup2 <- parApply(cl, grid,1, sup.rReg, N=N, me.response = me.response, ta = ta, tia = tia, tja = tja, tij = tij, T = T, topology = topology, times = times, model.type = model.type, ultrametric = ultrametric, Y = Y,  pred = pred, xx = xx, beta1 = beta1, error_condition = error_condition, s.X = s.X, n.pred = n.pred, num.prob = num.prob, cm2 = cm2, me.pred = me.pred, me.cov = me.cov, convergence = convergence, n.fixed = n.fixed)
        stopCluster(cl)
      } else {
        sup2 <- apply(grid, 1,sup.rReg, N=N, me.response = me.response, ta = ta, tia = tia, tja = tja, tij = tij, T = T, topology = topology, times = times, model.type = model.type, ultrametric = ultrametric, Y = Y,  pred = pred, xx = xx, beta1 = beta1, error_condition = error_condition, s.X = s.X, n.pred = n.pred, num.prob = num.prob, cm2 = cm2, me.pred = me.pred, me.cov = me.cov, convergence = convergence, n.fixed = n.fixed)
      }

      gof <- matrix(sup2, ncol=length(vy_values), byrow=TRUE, dimnames = list(half_life_values, vy_values))



      x<-rev(half_life_values)
      y<-vy_values
      z<-gof;

      ml<-max(z);
      for(i in 1:length(half_life_values))
      {
        for(j in 1:length(vy_values))
        {
          if(gof[i,j]==ml){alpha.est=log(2)/half_life_values[i]; vy.est=vy_values[j]}
        }
      }
      for(i in 1:length(half_life_values))
      {
        for(j in 1:length(vy_values))
        {
          if(gof[i,j]<=ml-support)gof[i, j]=ml-support;
        }
      }
      gof=gof-ml





      # FINAL OPTIMAL REGRESSION USING BEST ALPHA AND VY ESTIMATES #

      if(alpha.est==Inf)
      {
        gls.beta1<-glsyx.beta1<- solve(t(x.ols)%*%x.ols)%*%(t(x.ols)%*%Y)
        con.count<-0 # counter to break loop in the event of non-convergence



        repeat
        {
          s1<-as.numeric(s.X%*%(gls.beta1[2:(n.pred+1),]*gls.beta1[2:(n.pred+1),]))
          X<-cbind(1, pred)
          #V<-diag(rep(vy, times=N))+me.response+diag(as.numeric(me.pred%*%(gls.beta1[2:(n.pred+1),]*gls.beta1[2:(n.pred+1),])))-diag(as.numeric(me.cov%*%(2*gls.beta1[2:length(gls.beta1),])))

          #V<-diag(rep(vy, times=N))+na.exclude(me.response)+((Vu-(Vu%*%pseudoinverse(Vu+Vd)%*%Vu))*(gls.beta1[2:(n.pred+1),]*gls.beta1[2:(n.pred+1),]))-diag(as.numeric(me.cov%*%(2*gls.beta1[2:length(gls.beta1),])))



          V<-diag(rep(vy.est, times=N))+ na.exclude(me.response) + obs_var_con -diag(as.vector(na.exclude(me.cov%*%(2*gls.beta1[2:length(gls.beta1),]))));

          V.inverse<-solve(V)
          beta.i.var<-ev.beta.i.var<-pseudoinverse(t(X)%*%V.inverse%*%X)
          beta.i<-beta.i.var%*%(t(X)%*%V.inverse%*%Y)
          test<-matrix(nrow=(n.pred+1))
          for(f in 1:(n.pred+1))
          {
            if(abs(as.numeric(beta.i[f]-gls.beta1[f]))<=convergence) test[f]=0 else test[f]=1
          }
          if(sum(test)==0) break
          con.count=con.count+1
          if(con.count >= 50)
          {
            message("Warning, estimates did not converge after 50 iterations, last estimates printed out")
            break
          }
          gls.beta1<-glsyx.beta1<-beta.i
        }
        gls.beta1<-glsyx.beta1<-beta.i
        X<-cbind(1, pred)
        #V<-diag(rep(vy, times=N))+me.response+diag(as.numeric(me.pred%*%(gls.beta1[2:(n.pred+1),]*gls.beta1[2:(n.pred+1),])))-diag(as.numeric(me.cov%*%(2*gls.beta1[2:length(gls.beta1),])))



        V<-diag(rep(vy.est, times=N))+ na.exclude(me.response) + obs_var_con-diag(as.vector(na.exclude(me.cov%*%(2*gls.beta1[2:length(gls.beta1),]))));

        pred.mean<-X%*%gls.beta1
        g.mean<-(t(rep(1, times=N))%*%solve(V)%*%Y)/sum(solve(V));
        sst<-t(Y-g.mean)%*% solve(V)%*%(Y-g.mean)
        sse<-t(Y-pred.mean)%*%solve(V)%*%(Y-pred.mean)
        r.squared<-(sst-sse)/sst

        ###### Start of Bias correction ######
        adj<-matrix(data=0, ncol= 1+ n.pred, nrow=N)


        #if(ultrametric==TRUE)

        for(i in 1:n.pred)
        {
          adj[,(1+i)] <- as.numeric(sigma.X.estimate(pred[,i], me.pred[,i], topology, times)[1]);
        }


        if(ultrametric==TRUE) Vu<-Vu else Vu<-Vu[(length(Vu[1,])-(N*2)+1):length(Vu[1,]),(length(Vu[1,])-(N*2)+1):length(Vu[1,])]
        if(ultrametric==TRUE) Vd<-Vd else Vd<-Vd[(length(Vd[1,])-(N*2)+1):length(Vd[1,]),(length(Vd[1,])-(N*2)+1):length(Vd[1,])]



        correction<-matrix(Vu%*%pseudoinverse(Vd+Vu)%*%(c(X)-c(adj)),  ncol=ncol(X), nrow=nrow(X), byrow=F)



        bias_corr<-pseudoinverse(t(X)%*%V.inverse%*%X)%*%t(X)%*%V.inverse%*%correction

        m<-length(glsyx.beta1)



        corrected_betas<-solve(diag(1,m,m)-bias_corr)%*% glsyx.beta1

        ###### End of Bias correction ######
      }

      else
      {
        if(ultrametric==TRUE)
          gls.beta1<-solve(t(x.ols)%*%x.ols)%*%(t(x.ols)%*%Y)
        else
          gls.beta1<-rbind(0, 0, solve(t(x.ols)%*%x.ols)%*%(t(x.ols)%*%Y));
        con.count<-0;


        repeat
        {
          if(ultrametric==TRUE)
            s1<-as.numeric(s.X%*%(gls.beta1[2:(n.pred+1),]*gls.beta1[2:(n.pred+1),]))
          else
            s1<-as.numeric(s.X%*%(gls.beta1[4:(n.pred+3),]*gls.beta1[4:(n.pred+3),]))

          obs_var_con <-matrix(0, nrow=N, ncol=N)

          for (e in seq(from=1, to=ncol(x.ols), by=1)){
            for (j in seq(from=1, to=ncol(x.ols), by=1)) {
              tmp<-error_condition[xx[e]:(e*N),xx[j]:(j*N)]*(gls.beta1[e]*(1-(1-exp(-alpha.est*T))/(alpha.est*T)))*(gls.beta1[j]*(1-(1-exp(-alpha.est*T))/(alpha.est*T)))
              obs_var_con <-obs_var_con + tmp
            }

          }

          for(p in 1:N)
          {
            for(q in 1:N)
            {
              if(ta[q,p]==0)num.prob[q,p]=1 else num.prob[q,p]=(1-exp(-alpha.est*ta[q,p]))/(alpha.est*ta[q,p])
            }
          }
          cm1<-(s1/(2*alpha.est)+vy.est)*(1-exp(-2*alpha.est*ta))*exp(-alpha.est*tij)
          for(p in 1:N)
          {
            for(q in 1:N)
            {
              cm2[p,q]<-(((1-exp(-alpha.est*T[p]))/(alpha.est*T[p]))*((1-exp(-alpha.est*T[q]))/(alpha.est*T[q]))-(exp(-alpha.est*tia[p, q])*(1-exp(-alpha.est*T[p]))/(alpha.est*T[q])+ exp(-alpha.est*tja[p, q])*(1-exp(-alpha.est*T[p]))/(alpha.est*T[p]))*(num.prob[p,q]))
            }
          }
          if(ultrametric==TRUE)
          {
            X<-cbind(1, (1-(1-exp(-alpha.est*T))/(alpha.est*T))*pred)
            mv<-diag(rowSums(matrix(data=as.numeric(me.pred)*t(kronecker(gls.beta1[2:(n.pred+1), ], (1-(1-exp(-alpha.est*T))/(alpha.est*T)))^2), ncol=n.pred)))
            mcov<-diag(rowSums(matrix(data=as.numeric(me.cov)*t(kronecker(2*gls.beta1[2:length(gls.beta1),], (1-(1-exp(-alpha.est*T))/(alpha.est*T)))), ncol=n.pred)))

            #V<-cm1+(s1*ta*cm2)+me.response+mv-mcov;
            V<-cm1+(s1*ta*cm2)+na.exclude(me.response)+ obs_var_con-mcov

          }
          else
          {
            nu.X<-cbind(1-exp(-alpha.est*T), 1-exp(-alpha.est*T)-(1-(1-exp(-alpha.est*T))/(alpha.est*T)), exp(-alpha.est*T), (1-(1-exp(-alpha.est*T))/(alpha.est*T))*pred)
            mv<-diag(rowSums(matrix(data=as.numeric(me.pred)*t(kronecker(gls.beta1[4:(n.pred+3), ], (1-(1-exp(-alpha.est*T))/(alpha.est*T)))^2), ncol=n.pred)))

            mcov<-diag(rowSums(matrix(data=as.numeric(me.cov)*t(kronecker(2*gls.beta1[4:length(gls.beta1),], (1-(1-exp(-a*T))/(a*T)))), ncol=n.pred)));
            #V<-cm1+(s1*ta*cm2)+me.response+mv-mcov
            V<-cm1+(s1*ta*cm2)+na.exclude(me.response)+ obs_var_con-mcov

          }

          V.inverse<-solve(V)


          if(ultrametric==TRUE)
          {
            beta.i.var<-pseudoinverse(t(X)%*%V.inverse%*%X)
            beta.i<-beta.i.var%*%(t(X)%*%V.inverse%*%Y)
            test<-matrix(nrow=(n.pred+1))
            for(f in 1:(n.pred+1))
            {
              if(abs(as.numeric(beta.i[f]-gls.beta1[f]))<=convergence) test[f]=0 else test[f]=1
            }
            if(sum(test)==0) break
            con.count=con.count+1
            if(con.count >= 50)
            {
              message("Warning, estimates did not converge after 50 iterations, last estimates printed out")
              break
            }


            gls.beta1<-beta.i
          }
          else

          {

            beta.i.var<-pseudoinverse(t(nu.X)%*%V.inverse%*%nu.X)
            beta.i<-beta.i.var%*%(t(nu.X)%*%V.inverse%*%Y)
            test<-matrix(nrow=(n.pred))
            for(f in 4:(n.pred+3))
            {
              if(abs(as.numeric(beta.i[f]-gls.beta1[f]))<=convergence) test[(f-3)]=0 else test[(f-3)]=1
            }
            if(sum(test)==0) break
            con.count=con.count+1
            if(con.count >= 50)
            {
              message("Warning, estimates did not converge after 50 iterations, last estimates printed out")
              break
            }

            gls.beta1<-beta.i
          }
        }


        # END OF ITERATED GLS LOOP #


        # CODE FOR SST, SSE AND R-SQUARED #

        if(ultrametric==TRUE)
          gls.beta1<-beta.i
        else
        {
          gls.beta1<-beta.i
          ind.par<-matrix(data=0, nrow=N, ncol=4, dimnames=list(NULL, c("Bo", "Bi.Xia", "Yo", "Sum")))
          ind.par[,1]<-beta.i[1]*nu.X[,1]
          ind.par[,2]<-(beta.i[2]*nu.X[,2])
          ind.par[,3]<-beta.i[3]*nu.X[,3]
          ind.par[,4]<-ind.par[,1]+ind.par[,2]+ind.par[,3]
          mean.Bo=mean(ind.par[,4])
        }

        if(ultrametric==TRUE)
        {
          X<-cbind(1, (1-(1-exp(-alpha.est*T))/(alpha.est*T))*pred)

          obs_var_con <-matrix(0, nrow=N, ncol=N)

          for (e in seq(from=1, to=ncol(x.ols), by=1)){
            for (j in seq(from=1, to=ncol(x.ols), by=1)) {
              tmp<-error_condition[xx[e]:(e*N),xx[j]:(j*N)]*(gls.beta1[e]*(1-(1-exp(-alpha.est*T))/(alpha.est*T)))*(gls.beta1[j]*(1-(1-exp(-alpha.est*T))/(alpha.est*T)))
              obs_var_con <-obs_var_con + tmp
            }

          }


          mv<-diag(rowSums(matrix(data=as.numeric(me.pred)*t(kronecker(gls.beta1[2:(n.pred+1), ], (1-(1-exp(-alpha.est*T))/(alpha.est*T)))^2), ncol=n.pred)))
          mcov<-diag(rowSums(matrix(data=as.numeric(me.cov)*t(kronecker(2*gls.beta1[2:length(gls.beta1),], (1-(1-exp(-alpha.est*T))/(alpha.est*T)))), ncol=n.pred)))
          #V<-cm1+(s1*ta*cm2)+me.response+mv-mcov;
          V<-cm1+(s1*ta*cm2)+na.exclude(me.response)+ obs_var_con-mcov

          pred.mean<-X%*%gls.beta1

        }
        else
        {
          nu.X<-cbind(1-exp(-alpha.est*T), 1-exp(-alpha.est*T)-(1-(1-exp(-alpha.est*T))/(alpha.est*T)), exp(-alpha.est*T), (1-(1-exp(-alpha.est*T))/(alpha.est*T))*pred)

          obs_var_con <-matrix(0, nrow=N, ncol=N)

          for (e in seq(from=1, to=ncol(x.ols), by=1)){
            for (j in seq(from=1, to=ncol(x.ols), by=1)) {
              tmp<-error_condition[xx[e]:(e*N),xx[j]:(j*N)]*(gls.beta1[e]*(1-(1-exp(-alpha.est*T))/(alpha.est*T)))*(gls.beta1[j]*(1-(1-exp(-alpha.est*T))/(alpha.est*T)))
              obs_var_con <-obs_var_con + tmp
            }

          }


          mv<-diag(rowSums(matrix(data=as.numeric(me.pred)*t(kronecker(gls.beta1[4:(n.pred+3), ], (1-(1-exp(-alpha.est*T))/(alpha.est*T)))^2), ncol=n.pred)))
          mcov<-diag(rowSums(matrix(data=as.numeric(me.cov)*t(kronecker(2*gls.beta1[4:length(gls.beta1),], (1-(1-exp(-alpha.est*T))/(alpha.est*T)))), ncol=n.pred)))
          #V<-cm1+(s1*ta*cm2)+me.response+mv-mcov
          V<-cm1+(s1*ta*cm2)+na.exclude(me.response)+ obs_var_con-mcov


          pred.mean<-nu.X%*%gls.beta1
        }

        g.mean<-(t(rep(1, times=N))%*%solve(V)%*%Y)/sum(solve(V));
        sst<-t(Y-g.mean)%*% solve(V)%*%(Y-g.mean)
        sse<-t(Y-pred.mean)%*%solve(V)%*%(Y-pred.mean)
        r.squared<-(sst-sse)/sst


        # FINAL EVOLUTIONARY REGRESSION USING BEST ALPHA AND VY ESTIMATES AND KNOWN VARIANCE MATRIX #


        if(ultrametric==TRUE)  s1<-as.numeric(s.X%*%(gls.beta1[2:(n.pred+1),]*gls.beta1[2:(n.pred+1),]))
        else s1<-as.numeric(s.X%*%(gls.beta1[4:(n.pred+3),]*gls.beta1[4:(n.pred+3),]));

        obs_var_con <-matrix(0, nrow=N, ncol=N)

        for (e in seq(from=1, to=ncol(x.ols), by=1)){
          for (j in seq(from=1, to=ncol(x.ols), by=1)) {
            tmp<-error_condition[xx[e]:(e*N),xx[j]:(j*N)]*(gls.beta1[e]*(1-(1-exp(-alpha.est*T))/(alpha.est*T)))*(gls.beta1[j]*(1-(1-exp(-alpha.est*T))/(alpha.est*T)))
            obs_var_con <-obs_var_con + tmp
          }

        }


        for(p in 1:N)
        {
          for(q in 1:N)
          {
            if(ta[q,p]==0)num.prob[q,p]=1 else num.prob[q,p]=(1-exp(-alpha.est*ta[q,p]))/(alpha.est*ta[q,p])
          }
        }
        cm1<-(s1/(2*alpha.est)+vy.est)*(1-exp(-2*alpha.est*ta))*exp(-alpha.est*tij)
        for(p in 1:N)
        {
          for(q in 1:N)
          {
            cm2[p,q]<-(((1-exp(-alpha.est*T[p]))/(alpha.est*T[p]))*((1-exp(-alpha.est*T[q]))/(alpha.est*T[q]))-(exp(-alpha.est*tia[p, q])*(1-exp(-alpha.est*T[p]))/(alpha.est*T[q])+ exp(-alpha.est*tja[p, q])*(1-exp(-alpha.est*T[p]))/(alpha.est*T[p]))*(num.prob[p,q]))
          }
        }

        if(ultrametric==TRUE)

        {
          #V<-cm1+(s1*ta*cm2)+me.response+diag(as.numeric(me.pred%*%(gls.beta1[2:(n.pred+1),]*gls.beta1[2:(n.pred+1),])))-diag(as.numeric(me.cov%*%(2*gls.beta1[2:length(gls.beta1),])))


          V<-cm1+(s1*ta*cm2)+na.exclude(me.response)+ obs_var_con-diag(as.numeric(me.cov%*%(2*gls.beta1[2:length(gls.beta1),])))
        }


        else
          # V<-cm1+(s1*ta*cm2)+me.response+diag(as.numeric(me.pred%*%(gls.beta1[4:(n.pred+3),]*gls.beta1[4:(n.pred+3),])))-diag(as.numeric(me.cov%*%(2*gls.beta1[4:length(gls.beta1),])));

        {


          V<-cm1+(s1*ta*cm2)+na.exclude(me.response)+ obs_var_con-diag(as.numeric(me.cov%*%(2*gls.beta1[4:length(gls.beta1),])))
        }


        X1<-cbind(1, pred)
        V.inverse<-solve(V)
        ev.beta.i.var<-pseudoinverse(t(X1)%*%V.inverse%*%X1)
        ev.beta.i<-ev.beta.i.var%*%(t(X1)%*%V.inverse%*%Y)
        glsyx.beta1<-ev.beta.i


        ###### Start of Bias correction ######



        adj<-matrix(data=0, ncol= 1+ n.pred, nrow=N)


        #if(ultrametric==TRUE)

        for(i in 1:n.pred)
        {
          adj[,(1+i)] <- as.numeric(sigma.X.estimate(pred[,i], me.pred[,i], topology, times)[1]);
        }

        if(ultrametric==TRUE) Vu<-Vu else Vu<-Vu[(length(Vu[1,])-(N*2)+1):length(Vu[1,]),(length(Vu[1,])-(N*2)+1):length(Vu[1,])]
        if(ultrametric==TRUE) Vd<-Vd else Vd<-Vd[(length(Vd[1,])-(N*2)+1):length(Vd[1,]),(length(Vd[1,])-(N*2)+1):length(Vd[1,])]




        correction<-matrix(Vu%*%pseudoinverse(Vd+Vu)%*%(c(X1)-c(adj)),  ncol=ncol(X1), nrow=nrow(X1), byrow=F)



        bias_corr<-pseudoinverse(t(X1)%*%V.inverse%*%X1)%*%t(X1)%*%V.inverse%*%correction

        m<-length(glsyx.beta1)



        corrected_betas<-solve(diag(1,m,m)-bias_corr)%*% glsyx.beta1


        ###### End of Bias correction ######








      }
      # END OF HALFLIFE 0 CONDITION #



    } # END OF RANDOM COVARIATE REGRESSION ESTIMATION

  }# END OF FIXED COVARIATE, MIXED OR RANDOM MODELS PARAMETER ESTIMATION



  # EVALUATE IF IT IS A FIXED AND RANDOM COVARIATE ANCOVA OR REGRESSION MODEL ESTIMATE PARAMETERS WITH ITERATED GLS TO A) TAKE MEASUREMENT VARIANCE INTO ACCOUNT OR B) RANDOM EFFECTS INTO ACCOUNT IN THE CASE OF THE MIXED MODEL AND REGRESSION

  if(model.type == "mmfANCOVA" || model.type=="mfReg")
  {
    # SET UP INITIAL MATRICES FOR MULTIPLE REGRESSION AND CALCULATE THETA AND SIGMA FOR RANDOM PREDICTOR / S

    pred<-data.frame(random.cov);
    n.pred<-length(pred[1,]);
    pred<-matrix(data=pred[!is.na(pred)], ncol=n.pred);
    if(is.null(me.random.cov)) me.pred<-matrix(data=0, nrow=N, ncol=n.pred) else me.pred<-matrix(data=me.random.cov[!is.na(me.random.cov)], ncol=n.pred);
    if(is.null(mecov.random.cov)) me.cov<-matrix(data=0, nrow=N, ncol=n.pred) else me.cov<-matrix(data=mecov.random.cov[!is.na(mecov.random.cov)], ncol=n.pred);

    s.X<-matrix(data=0, ncol=n.pred)  # PREDICTOR SIGMA
    for(i in 1:n.pred)
    {
      s.X[,i] <- as.numeric(sigma.X.estimate(pred[,i],me.pred[,i], topology, times)[2]);
    }

    theta.X<-matrix(data=0, ncol=n.pred)  #PREDICTOR THETA
    for(i in 1:n.pred)
    {
      theta.X[,i] <- as.numeric(sigma.X.estimate(pred[,i],me.pred[,i], topology, times)[1]);
    }



    # END OF RANDOM PREDICTOR THETA AND SIGMA ESTIMATES

    # FIXED COVARIATES

    fixed.pred<-data.frame(fixed.cov);
    n.fixed.pred<-length(fixed.pred[1,]);
    fixed.pred<-matrix(data=fixed.pred[!is.na(fixed.pred)], ncol=n.fixed.pred);
    if(is.null(me.fixed.cov)) me.fixed.pred<-matrix(data=0, nrow=N, ncol=n.fixed.pred) else me.fixed.pred<- matrix(data=me.fixed.cov[!is.na(me.fixed.cov)], ncol=n.fixed.pred);
    if(is.null(mecov.fixed.cov)) me.fixed.cov<-matrix(data=0, nrow=N, ncol=n.fixed.pred) else me.fixed.cov<-matrix(data=me.cov.fixed.cov[!is.na(me.cov.fixed.cov)], ncol=n.fixed.pred);


    ## INITIAL OLS ESTIMATES TO SEED ITERATED GLS

    if(model.type=="mfReg")
    {
      x.ols<-cbind(1, fixed.pred, pred);
      beta1<-solve(t(x.ols)%*%x.ols)%*%(t(x.ols)%*%Y);
      if(ultrametric == FALSE) beta1<-rbind(0, 0, beta1); # 2 additional parameter seeds for Ya and Xa


      #Defining the dimensionality of Vd
      Vd<-matrix(0,ncol=(N*length(beta1[,1])), nrow=(N*length(beta1[,1])))


      #Putting in elements in VD for fixed covariates
      true_var<-matrix(data=0, ncol=n.fixed.pred, nrow=N);
      for (i in 1:n.fixed.pred)
      {
        true_var[,i]<-var(na.exclude(fixed.pred[,i]))-as.numeric(na.exclude(me.fixed.pred[,i]))
      }

      true_var<-c(true_var)
      if(ultrametric == TRUE) Vd[1:((N*n.fixed.pred)+N),1:((N*n.fixed.pred)+N)]<-diag(c(rep(0,N),true_var)) else Vd[((3*N+1):((3*N+(n.fixed.pred*N)))),((3*N+1):((3*N+(n.fixed.pred*N))))]<-diag(c(true_var))


      #Putting in elements in VD for random covariates

      xx<-seq(from=1, to	=length(Vd[,1]), by=N)

      if(ultrametric == TRUE) xx<-xx[-(1:(1+ n.fixed.pred))] else xx<-xx[-(1:(3+ n.fixed.pred))]

      yy<-seq(from=N, to	=length(Vd[,1]), by=N)

      if(ultrametric == TRUE) yy<-yy[-(1:(1+ n.fixed.pred))] else yy<-yy[-(1:(3+ n.fixed.pred))]


      for (i in seq(from=1, to=nrow(s.X), by=1)){
        Vd[xx[i]:yy[i],xx[i]:yy[i]]<-pt$bt*s.X[,i]

      }


      # Defining Vu
      if(ultrametric == TRUE) Vu<-diag(c(rep(0,N), c(as.numeric(na.exclude(me.fixed.pred))),c(as.numeric(na.exclude(me.pred))))) else Vu<-diag(c(rep(0,N*3), c(as.numeric(na.exclude(me.pred))), c(as.numeric(na.exclude(me.fixed.pred)))))

      error_condition<-Vu-(Vu%*%pseudoinverse(Vu+Vd)%*%Vu)

      mybiglist <- list()

      xx<-seq(from=1, to=length(Vu[,1]), by=N)
      mybiglist <- list()


      for (i in seq(from=1, to=nrow(beta1), by=1)){
        for (j in seq(from=1, to=nrow(beta1), by=1)) {
          tmp <- list(error_condition[xx[i]:(i*N),xx[j]:(j*N)]*beta1[i]*beta1[j])
          mybiglist[xx[i]+j] <- tmp
        }

      }


      mybiglist <-rmNullObs(mybiglist)
      obs_var_con<-Reduce('+', mybiglist)


    }

    if(model.type=="mmfANCOVA")
    {
      regime.specs<-fixed.fact;
      n.fixed<-length(levels(as.factor(regime.specs)))
      regime.specs<-as.factor(regime.specs)
      x.ols<-cbind(weight.matrix(1000000000000000000000, topology, times, N, regime.specs, fixed.pred, intercept), pred);
      beta1<-solve(t(x.ols)%*%x.ols)%*%(t(x.ols)%*%Y);



      #Defining the dimensionality of Vd ## Might need to change the dimensionality of this matrix for different versions of the model and how that impacts the number of columns in x.ols
      Vd<-matrix(0,ncol=(N*length(beta1[,1])), nrow=(N*length(beta1[,1])))


      #Putting in elements in VD for fixed covariates
      true_var<-matrix(data=0, ncol=n.fixed.pred, nrow=N);
      for (i in 1:n.fixed.pred)
      {
        true_var[,i]<-var(na.exclude(fixed.pred[,i]))-as.numeric(na.exclude(me.fixed.pred[,i]))
      }

      true_var<-c(true_var)
      Vd[(((N* n.fixed)+(1)):(((N* n.fixed)+(1))+((N*n.fixed.pred)-1))),(((N* n.fixed)+(1)):(((N* n.fixed)+(1))+((N*n.fixed.pred)-1)))]<-diag(c(true_var))


      #Putting in elements in VD for random covariates

      xx<-seq(from=1, to	=length(Vd[,1]), by=N)

      if(ultrametric == TRUE) xx<-xx[-(1:(n.fixed+n.fixed.pred))] else xx<-xx[-(1:(n.fixed+ n.fixed.pred))]

      yy<-seq(from=N, to	=length(Vd[,1]), by=N)

      if(ultrametric == TRUE) yy<-yy[-(1:(n.fixed+n.fixed.pred))] else yy<-yy[-(1:(n.fixed+n.fixed.pred))]


      for (i in seq(from=1, to=nrow(s.X), by=1)){
        Vd[xx[i]:yy[i],xx[i]:yy[i]]<-pt$bt*s.X[,i]

      }

      # Defining Vu
      if(ultrametric == TRUE) Vu<-diag(c(rep(0,(N*(n.fixed))), c(as.numeric(na.exclude(me.fixed.pred))),c(as.numeric(na.exclude(me.pred))))) else Vu<-diag(c(rep(0,N*(2+ n.fixed))), c(as.numeric(na.exclude(me.fixed.pred))), c(as.numeric(na.exclude(me.pred))))

      error_condition<-Vu-(Vu%*%pseudoinverse(Vu+Vd)%*%Vu)

      mybiglist <- list()

      xx<-seq(from=1, to=length(Vu[,1]), by=N)
      mybiglist <- list()


      for (i in seq(from=1, to=nrow(beta1), by=1)){
        for (j in seq(from=1, to=nrow(beta1), by=1)) {
          tmp <- list(error_condition[xx[i]:(i*N),xx[j]:(j*N)]*beta1[i]*beta1[j])
          mybiglist[xx[i]+j] <- tmp
        }

      }


      mybiglist <-rmNullObs(mybiglist)
      obs_var_con<-Reduce('+', mybiglist)


    }



    # GRID ESTIMATION ROUTINE AND ITERATED GLS FOR MODELS THAT INCLUDE RANDOM EFFECTS

    if(model.type=="mmfANCOVA")
    {

      cat(c("   ", "t1/2  ", "Vy    ", "Supp  ", GS_head, if(is.null(dim(fixed.cov))) deparse(substitute(fixed.cov)) else colnames(fixed.cov), if(is.null(dim(random.cov))) deparse(substitute(random.cov)) else colnames(random.cov)), sep="   ");


      message(" ");

      for(i in 1:length(half_life_values))
      {
        for(k in 1:length(vy_values))
        {
          if(half_life_values[i]==0) a<-1000000000000000000000 else a <- ln2/half_life_values[i];
          vy <- vy_values[k];
          X<-cbind(weight.matrix(a, topology, times, N, regime.specs, fixed.pred, intercept), (1-(1-exp(-a*T))/(a*T))*pred);
          if(length(X[1,]) > length(beta1)) {beta1<-as.matrix(c(0, beta1)); n.fixed<-n.fixed+1}
          if(length(X[1,])< length(beta1)) {beta1<-solve(t(x.ols)%*%x.ols)%*%(t(x.ols)%*%Y);n.fixed<-length(levels(as.factor(regime.specs))); print("The Ya parameter is dropped as its coefficient is too small");}

          # CODE FOR ESTIMATING BETA USING ITERATED GLS

          con.count<-0;  # Counter for loop break if Beta's dont converge #
          repeat
          {
            if(half_life_values[i]==0)
            {
              X<-cbind(weight.matrix(1000000000000000000000, topology, times, N, regime.specs, fixed.pred, intercept), pred);
              #V<-diag(rep(vy, times=N))+me.response+diag(as.numeric(me.pred%*%(beta1[(n.fixed+1+n.fixed.pred):length(beta1),]*beta1[(n.fixed+1+n.fixed.pred):length(beta1),])))-diag(as.numeric(me.cov%*%(2*beta1[(n.fixed+1+n.fixed.pred):length(beta1),]))) + diag(as.numeric(me.fixed.pred%*%(beta1[(n.fixed+1):(length(beta1)-n.pred),]*beta1[(n.fixed+1):(length(beta1)-n.pred),])))-diag(as.numeric(me.fixed.cov%*%(2*beta1[(n.fixed+1):(length(beta1)-n.pred),])));

              V<-diag(rep(vy, times=N))+na.exclude(me.response)+ obs_var_con-diag(as.numeric(me.cov%*%(2*beta1[(n.fixed+1+n.fixed.pred):length(beta1),])))-diag(as.numeric(me.fixed.cov%*%(2*beta1[(n.fixed+1):(length(beta1)-n.pred),])));

            }
            else
            {

              X<-cbind(weight.matrix(a, topology, times, N, regime.specs, fixed.pred, intercept), (1-(1-exp(-a*T))/(a*T))*pred);


              s1<-as.numeric(s.X%*%(beta1[(n.fixed+1+n.fixed.pred):length(beta1),]*beta1[(n.fixed+1+n.fixed.pred):length(beta1),]));


              for(p in 1:N)
              {
                for(q in 1:N)
                {
                  if(ta[q,p]==0)num.prob[q,p]=1 else num.prob[q,p]=(1-exp(-a*ta[q,p]))/(a*ta[q,p]);
                }
              }
              cm1<-(s1/(2*a)+vy)*(1-exp(-2*a*ta))*exp(-a*tij);
              for(p in 1:N)
              {
                for(q in 1:N)
                {
                  cm2[p,q]<-(((1-exp(-a*T[p]))/(a*T[p]))*((1-exp(-a*T[q]))/(a*T[q]))-(exp(-a*tia[p, q])*(1-exp(-a*T[p]))/ (a*T[q])+ exp(-a*tja[p, q])*(1-exp(-a*T[p]))/(a*T[p]))*(num.prob[p,q]));
                }
              }

              mv<-diag(rowSums(matrix(data=as.numeric(me.pred)*t(kronecker(beta1[(n.fixed+1+n.fixed.pred):length(beta1), ], (1-(1-exp(-a*T))/(a*T)))^2), ncol=n.pred)));
              mcov<-diag(rowSums(matrix(data=as.numeric(me.cov)*t(kronecker(2*beta1[(n.fixed+1+n.fixed.pred):length(beta1),], (1-(1-exp(-a*T))/(a*T)))), ncol=n.pred)));
              mv.fixed<-diag(rowSums(matrix(data=as.numeric(me.fixed.pred)*t(kronecker(beta1[(n.fixed+1):(length(beta1)-n.pred), ], rep(1, times=N))), ncol=n.fixed.pred)));
              mcov.fixed<-diag(rowSums(matrix(data=as.numeric(me.fixed.cov)*t(kronecker(2*beta1[(n.fixed+1):(length(beta1)-n.pred),], rep(1, times=N))), ncol=n.fixed.pred)));



              #V<-cm1+(s1*ta*cm2)+me.response+mv+ mv.fixed-mcov-mcov.fixed;

              V<-cm1+(s1*ta*cm2)+na.exclude(me.response)+obs_var_con-mcov-mcov.fixed;


              obs_var_con
            } # END OF If ELSE CONDITION FOR HALF-LIFE 0 OR NOT

            # INTERMEDIATE ESTIMATION OF OPTIMAL REGRESSION #

            V.inverse<-solve(V)
            if(half_life_values[i]==0)
            {
              beta.i<-pseudoinverse(t(X)%*%V.inverse%*%X)%*%(t(X)%*%V.inverse%*%Y)
              test<-matrix(nrow=(length(beta.i)))
              for(f in 1:(length(beta.i)))
              {
                if(abs(as.numeric(beta.i[f]-beta1[f]))<=convergence) test[f]=0 else test[f]=1
              }
              if(sum(test)==0) break
              con.count=con.count+1
              if(con.count >= 50)
              {
                message("Warning, estimates did not converge after 50 iterations, last estimates printed out")
                break
              }

              beta1<-beta.i
            }
            else
            {

              beta.i<-pseudoinverse(t(X)%*%V.inverse%*%X)%*%(t(X)%*%V.inverse%*%Y)
              test<-matrix(nrow=(length(beta.i)))
              for(f in 1:(length(beta.i)))
              {
                if(abs(as.numeric(beta.i[f]-beta1[f]))<=convergence) test[f]=0 else test[f]=1
              }
              if(sum(test)==0) break
              con.count=con.count+1
              if(con.count >= 50)
              {
                message("Warning, estimates did not converge after 50 iterations, last estimates printed out")
                break
              }

              beta1<-beta.i
            }
          }


          ### END OF ITERATED GLS ESTIMATION FOR BETA #

          if(half_life_values[i]==0)
          {
            X<-cbind(weight.matrix(1000000000000000000000, topology, times, N, regime.specs, fixed.pred, intercept), pred)
            #V<-diag(rep(vy, times=N))+me.response+diag(as.numeric(me.pred%*%(beta1[(n.fixed+1+n.fixed.pred):length(beta1),]*beta1[(n.fixed+1+n.fixed.pred):length(beta1),])))-diag(as.numeric(me.cov%*%(2*beta1[(n.fixed+1+n.fixed.pred):length(beta1),]))) + diag(as.numeric(me.fixed.pred%*%(beta1[(n.fixed+1):(length(beta1)-n.pred),]*beta1[(n.fixed+1):(length(beta1)-n.pred),])))-diag(as.numeric(me.fixed.cov%*%(2*beta1[(n.fixed+1):(length(beta1)-n.pred),])));

            V<-diag(rep(vy, times=N))+na.exclude(me.response)+ obs_var_con-diag(as.numeric(me.cov%*%(2*beta1[(n.fixed+1+n.fixed.pred):length(beta1),])))-diag(as.numeric(me.fixed.cov%*%(2*beta1[(n.fixed+1):(length(beta1)-n.pred),])));


            V.inverse<-solve(V)
            eY<-X%*%beta1
            resid<-Y-eY;
            gof[i, k] <- -N/2*log(2*pi)-0.5*log(det(V))-0.5*(t(resid) %*% V.inverse%*%resid);
          }
          else
          {
            s1<-as.numeric(s.X%*%(beta1[(n.fixed+1+n.fixed.pred):length(beta1),]*beta1[(n.fixed+1+n.fixed.pred):length(beta1),]));
            for(p in 1:N)
            {
              for(q in 1:N)
              {
                if(ta[q,p]==0)num.prob[q,p]=1 else num.prob[q,p]=(1-exp(-a*ta[q,p]))/(a*ta[q,p]);
              }
            }
            cm1<-(s1/(2*a)+vy)*(1-exp(-2*a*ta))*exp(-a*tij);
            for(p in 1:N)
            {
              for(q in 1:N)
              {
                cm2[p,q]<-(((1-exp(-a*T[p]))/(a*T[p]))*((1-exp(-a*T[q]))/(a*T[q]))-(exp(-a*tia[p, q])*(1-exp(-a*T[p]))/(a*T[q])+ exp(-a*tja[p, q])*(1-exp(-a*T[p]))/(a*T[p]))*(num.prob[p,q]));
              }
            }
            X<-cbind(weight.matrix(a, topology, times, N, regime.specs, fixed.pred, intercept), (1-(1-exp(-a*T))/(a*T))*pred);
            mv<-diag(rowSums(matrix(data=as.numeric(me.pred)*t(kronecker(beta1[(n.fixed+1+n.fixed.pred):length(beta1), ], (1-(1-exp(-a*T))/(a*T)))^2), ncol=n.pred)));
            mcov<-diag(rowSums(matrix(data=as.numeric(me.cov)*t(kronecker(2*beta1[(n.fixed+1+n.fixed.pred):length(beta1),], (1-(1-exp(-a*T))/(a*T)))), ncol=n.pred)));
            mv.fixed<-diag(rowSums(matrix(data=as.numeric(me.fixed.pred)*t(kronecker(beta1[(n.fixed+1):(length(beta1)-n.pred), ], rep(1, times=N))), ncol=n.fixed.pred)));
            mcov.fixed<-diag(rowSums(matrix(data=as.numeric(me.fixed.cov)*t(kronecker(2*beta1[(n.fixed+1):(length(beta1)-n.pred),], rep(1, times=N))), ncol=n.fixed.pred)));

            V<-cm1+(s1*ta*cm2)+na.exclude(me.response)+ obs_var_con -mcov-mcov.fixed;
            V.inverse<-solve(V)

            eY<-X%*%beta1

            resid<-Y-eY;
            gof[i, k] <- -N/2*log(2*pi)-0.5*log(det(V))-0.5*(t(resid) %*% V.inverse%*%resid);
          }  # END OF CONDITION FOR HALF-LIFE = 0 #
          print(as.numeric(round(cbind(if(a!=0)log(2)/a else 0.00, vy, gof[i,k], t(beta1)), 4)))
        }
      }



      # END OF GRID SETUP,START OF GRID SEARCH FOR BEST ALPHA AND VY ESTIMATES #

      x<-rev(half_life_values)
      y<-vy_values
      z<-gof;
      ml<-max(z);
      for(i in 1:length(half_life_values))
      {
        for(j in 1:length(vy_values))
        {
          if(gof[i,j]==ml){alpha.est=log(2)/half_life_values[i]; vy.est=vy_values[j]}
        }
      }
      for(i in 1:length(half_life_values))
      {
        for(j in 1:length(vy_values))
        {
          if(gof[i,j]<=ml-support)gof[i, j]=ml-support;
        }
      }
      gof=gof-ml


      n.fixed<-length(levels(as.factor(regime.specs)))   ### reset before final regression


      # FINAL OPTIMAL REGRESSION USING BEST ALPHA AND VY ESTIMATES #

      if(alpha.est==Inf || alpha.est >=1000000000000000000000)
      {
        x.ols<-cbind(weight.matrix(1000000000000000000000, topology, times, N, regime.specs, fixed.pred, intercept), pred)
        gls.beta1<-solve(t(x.ols)%*%x.ols)%*%(t(x.ols)%*%Y)
        con.count<-0;
        repeat
        {

          s1<-as.numeric(s.X%*%(gls.beta1[(n.fixed+1+n.fixed.pred):length(gls.beta1),]*gls.beta1[(n.fixed+1+n.fixed.pred):length(gls.beta1),]));
          X<-cbind(weight.matrix(1000000000000000000000, topology, times, N, regime.specs, fixed.cov, intercept), pred)
          #V<-diag(rep(vy, times=N))+me.response+diag(as.numeric(me.pred%*%(gls.beta1[(n.fixed+1+n.fixed.pred):length(gls.beta1),]*gls.beta1[(n.fixed+1+n.fixed.pred):length(gls.beta1),])))-diag(as.numeric(me.cov%*%(2*gls.beta1[(n.fixed+1+n.fixed.pred):length(gls.beta1),]))) + diag(as.numeric(me.fixed.pred%*%(gls.beta1[(n.fixed+1):(length(gls.beta1)-n.pred),]*gls.beta1[(n.fixed+1):(length(gls.beta1)-n.pred),])))-diag(as.numeric(me.fixed.cov%*%(2*gls.beta1[(n.fixed+1):(length(gls.beta1)-n.pred),])));

          V<-diag(rep(vy, times=N))+na.exclude(me.response)+ obs_var_con -diag(as.numeric(me.cov%*%(2*gls.beta1[(n.fixed+1+n.fixed.pred):length(gls.beta1),]))) -diag(as.numeric(me.fixed.cov%*%(2*gls.beta1[(n.fixed+1):(length(gls.beta1)-n.pred),])));

          V.inverse<-solve(V)
          beta.i.var<-ev.beta.i.var<-pseudoinverse(t(X)%*%V.inverse%*%X)
          beta.i<-beta.i.var%*%(t(X)%*%V.inverse%*%Y)
          test<-matrix(nrow=(length(beta.i)))
          for(f in 1:(length(beta.i)))
          {
            if(abs(as.numeric(beta.i[f]-gls.beta1[f]))<=convergence) test[f]=0 else test[f]=1
          }
          if(sum(test)==0) break
          con.count=con.count+1
          if(con.count >= 50)
          {
            message("Warning, estimates did not converge after 50 iterations, last estimates printed out")
            break
          }
          gls.beta1<-beta.i
        }
        gls.beta1<-beta.i
        X<-cbind(weight.matrix(1000000000000000000000, topology, times, N, regime.specs, fixed.pred, intercept), pred)
        #V<-diag(rep(vy, times=N))+me.response+diag(as.numeric(me.pred%*%(gls.beta1[(n.fixed+1+n.fixed.pred):length(gls.beta1),]*gls.beta1[(n.fixed+1+n.fixed.pred):length(gls.beta1),])))-diag(as.numeric(me.cov%*%(2*gls.beta1[(n.fixed+1+n.fixed.pred):length(gls.beta1),]))) + diag(as.numeric(me.fixed.pred%*%(gls.beta1[(n.fixed+1):(length(gls.beta1)-n.pred),]*gls.beta1[(n.fixed+1):(length(gls.beta1)-n.pred),])))-diag(as.numeric(me.fixed.cov%*%(2*gls.beta1[(n.fixed+1):(length(gls.beta1)-n.pred),])));

        V<-diag(rep(vy, times=N))+na.exclude(me.response)+  obs_var_con -diag(as.numeric(me.cov%*%(2*gls.beta1[(n.fixed+1+n.fixed.pred):length(gls.beta1),]))) -diag(as.numeric(me.fixed.cov%*%(2*gls.beta1[(n.fixed+1):(length(gls.beta1)-n.pred),])));

        pred.mean<-X%*%gls.beta1
        g.mean<-(t(rep(1, times=N))%*%solve(V)%*%Y)/sum(solve(V));
        sst<-t(Y-g.mean)%*% solve(V)%*%(Y-g.mean)

        sse<-t(Y-pred.mean)%*%solve(V)%*%(Y-pred.mean)

        r.squared<-(sst-sse)/sst


      }
      else
      {
        x.ols<-cbind(weight.matrix(1000000000000000000000, topology, times, N, regime.specs, fixed.pred, intercept), pred)
        gls.beta1<-solve(t(x.ols)%*%x.ols)%*%(t(x.ols)%*%Y)
        con.count<-0;

        X<-cbind(weight.matrix(alpha.est, topology, times, N, regime.specs, fixed.pred, intercept), (1-(1-exp(-alpha.est*T))/(alpha.est*T))*pred);
        if(length(X[1,]) > length(gls.beta1)) {gls.beta1<-as.matrix(c(0, gls.beta1)); n.fixed<-n.fixed+1}
        if(length(X[1,])< length(gls.beta1)) {gls.beta1<-solve(t(x.ols)%*%x.ols)%*%(t(x.ols)%*%Y);n.fixed<-length(levels(as.factor(regime.specs)))}
        repeat
        {

          X<-cbind(weight.matrix(alpha.est, topology, times, N, regime.specs, fixed.pred, intercept), (1-(1-exp(-alpha.est*T))/(alpha.est*T))*pred);
          s1<-as.numeric(s.X%*%(gls.beta1[(n.fixed+1+n.fixed.pred):length(gls.beta1),]*gls.beta1[(n.fixed+1+n.fixed.pred):length(gls.beta1),]));


          for(p in 1:N)
          {
            for(q in 1:N)
            {
              if(ta[q,p]==0)num.prob[q,p]=1 else num.prob[q,p]=(1-exp(-alpha.est*ta[q,p]))/(alpha.est*ta[q,p])
            }
          }
          cm1<-(s1/(2*alpha.est)+vy.est)*(1-exp(-2*alpha.est*ta))*exp(-alpha.est*tij)
          for(p in 1:N)
          {
            for(q in 1:N)
            {
              cm2[p,q]<-(((1-exp(-alpha.est*T[p]))/(alpha.est*T[p]))*((1-exp(-alpha.est*T[q]))/(alpha.est*T[q]))-(exp(-alpha.est*tia[p, q])*(1-exp(-alpha.est*T[p]))/(alpha.est*T[q])+ exp(-alpha.est*tja[p, q])*(1-exp(-alpha.est*T[p]))/(alpha.est*T[p]))*(num.prob[p,q]))
            }
          }


          mv<-diag(rowSums(matrix(data=as.numeric(me.pred)*t(kronecker(gls.beta1[(n.fixed+1+n.fixed.pred):length(gls.beta1), ], (1-(1-exp(-alpha.est*T))/(alpha.est*T)))^2), ncol=n.pred)));
          mcov<-diag(rowSums(matrix(data=as.numeric(me.cov)*t(kronecker(2*gls.beta1[(n.fixed+1+n.fixed.pred):length(gls.beta1),], (1-(1-exp(-alpha.est*T))/(alpha.est*T)))), ncol=n.pred)));
          mv.fixed<-diag(rowSums(matrix(data=as.numeric(me.fixed.pred)*t(kronecker(gls.beta1[(n.fixed+1):(length(gls.beta1)-n.pred), ], rep(1, times=N))), ncol=n.fixed.pred)));
          mcov.fixed<-diag(rowSums(matrix(data=as.numeric(me.fixed.cov)*t(kronecker(2*gls.beta1[(n.fixed+1):(length(gls.beta1)-n.pred),], rep(1, times=N))), ncol=n.fixed.pred)));

          #V<-cm1+(s1*ta*cm2)+me.response+mv+ mv.fixed-mcov-mcov.fixed;
          V<-cm1+(s1*ta*cm2)+na.exclude(me.response)+ obs_var_con-mcov-mcov.fixed;


          V.inverse<-solve(V)
          beta.i.var<-pseudoinverse(t(X)%*%V.inverse%*%X)
          beta.i<-beta.i.var%*%(t(X)%*%V.inverse%*%Y)
          test<-matrix(nrow=(length(beta.i)))
          for(f in 1:(length(beta.i)))
          {
            if(abs(as.numeric(beta.i[f]-gls.beta1[f]))<=convergence) test[f]=0 else test[f]=1
          }
          if(sum(test)==0) break
          con.count=con.count+1
          if(con.count >= 50)
          {
            message("Warning, estimates did not converge after 50 iterations, last estimates printed out")
            break
          }


          gls.beta1<-beta.i

          mybiglist <- list()

          xx<-seq(from=1, to=length(Vu[,1]), by=N)
          mybiglist <- list()


          for (i in seq(from=1, to=nrow(gls.beta1), by=1)){
            for (j in seq(from=1, to=nrow(gls.beta1), by=1)) {
              tmp <- list(error_condition[xx[i]:(i*N),xx[j]:(j*N)]* gls.beta1[i]* gls.beta1[j])
              mybiglist[xx[i]+j] <- tmp
            }

          }


          mybiglist <-rmNullObs(mybiglist)
          obs_var_con<-Reduce('+', mybiglist)


          X<-cbind(weight.matrix(alpha.est, topology, times, N, regime.specs, fixed.pred, intercept), (1-(1-exp(-alpha.est*T))/(alpha.est*T))*pred)

          mv<-diag(rowSums(matrix(data=as.numeric(me.pred)*t(kronecker(gls.beta1[(n.fixed+1+n.fixed.pred):length(gls.beta1), ], (1-(1-exp(-alpha.est*T))/(alpha.est*T)))^2), ncol=n.pred)));
          mcov<-diag(rowSums(matrix(data=as.numeric(me.cov)*t(kronecker(2*gls.beta1[(n.fixed+1+n.fixed.pred):length(gls.beta1),], (1-(1-exp(-alpha.est*T))/(alpha.est*T)))), ncol=n.pred)));
          mv.fixed<-diag(rowSums(matrix(data=as.numeric(me.fixed.pred)*t(kronecker(gls.beta1[(n.fixed+1):(length(gls.beta1)-n.pred), ], rep(1, times=N))), ncol=n.fixed.pred)));
          mcov.fixed<-diag(rowSums(matrix(data=as.numeric(me.fixed.cov)*t(kronecker(2*gls.beta1[(n.fixed+1):(length(gls.beta1)-n.pred),], rep(1, times=N))), ncol=n.fixed.pred)));

          #V<-cm1+(s1*ta*cm2)+me.response+mv+ mv.fixed-mcov-mcov.fixed;
          V<-cm1+(s1*ta*cm2)+na.exclude(me.response) + obs_var_con-mcov-mcov.fixed;


          pred.mean<-X%*%gls.beta1
          g.mean<-(t(rep(1, times=N))%*%solve(V)%*%Y)/sum(solve(V));
          sst<-t(Y-g.mean)%*% solve(V)%*%(Y-g.mean)
          sse<-t(Y-pred.mean)%*%solve(V)%*%(Y-pred.mean)
          r.squared<-(sst-sse)/sst



        }
      }

      # END OF ITERATED GLS LOOP #

      ###### Start of Bias correction ######



      adj<-matrix(data=0, ncol= ncol(X), nrow=N)

      for(i in 1:(n.fixed.pred))
      {
        adj[,(n.fixed+i)]<-mean(X[,(n.fixed+i)])
      }

      for(i in 1:(n.pred))
      {
        adj[,(n.fixed+ n.fixed.pred +i)] <- as.numeric(sigma.X.estimate(pred[,i], me.pred[,i], topology, times)[1]);
      }

      V.inverse<-solve(V)

      correction<-matrix(Vu%*%pseudoinverse(Vd+Vu)%*%(c(X)-c(adj)),  ncol=ncol(X), nrow=nrow(X), byrow=F)



      bias_corr<-pseudoinverse(t(X)%*%V.inverse%*%X)%*%t(X)%*%V.inverse%*%correction


      m<-length(gls.beta1);



      K<-solve(diag(1,m,m)-bias_corr)


      corrected_betas<-solve(diag(1,m,m)-bias_corr)%*% gls.beta1



      ###### End of Bias correction ######






    } # END OF ESTIMATION MIXED MODEL ANCOVA



    if(model.type=="mfReg")
    {

      if(ultrametric==TRUE)
      {

        cat(c("   ", "t1/2  ", "Vy    ", "Supp  ", "K     ",if(is.null(dim(fixed.cov))) deparse(substitute(fixed.cov)) else colnames(fixed.cov), if(is.null(dim(random.cov))) deparse(substitute(random.cov)) else colnames(random.cov)), sep="   ");

      }

      else
        cat(c("   ", "t1/2  ", "Vy    ", "Supp  ", "Ya    ", "Xa    " ,"Bo    ", if(is.null(dim(fixed.cov))) deparse(substitute(fixed.cov)) else colnames(fixed.cov), if(is.null(dim(random.cov))) deparse(substitute(random.cov)) else colnames(random.cov)), sep="   ");

      message(" ");

      for(i in 1:length(half_life_values))
      {
        for(k in 1:length(vy_values))
        {
          if(half_life_values[i]==0)
          {
            x.ols<-cbind(1, fixed.pred, pred)
            beta1<-solve(t(x.ols)%*%x.ols)%*%(t(x.ols)%*%Y)
            vy <- vy_values[k];
          }
          else
          {
            a <- ln2/half_life_values[i];
            vy <- vy_values[k];
            x.ols<-cbind(1,fixed.pred, pred)
            if(ultrametric==TRUE)
              beta1<-solve(t(x.ols)%*%x.ols)%*%(t(x.ols)%*%Y)
            else
              beta1<-rbind(0, 0, solve(t(x.ols)%*%x.ols)%*%(t(x.ols)%*%Y))
          }

          ### CODE FOR ESTIMATING BETA USING ITERATED GLS ###
          con.count<-0;  # Counter for loop break if Beta's dont converge #
          repeat
          {
            if(half_life_values[i]==0)
            {
              a<-Inf
              s1<-as.numeric(s.X%*%(beta1[(2+n.fixed.pred):(n.pred+1+n.fixed.pred),]*beta1[(2+n.fixed.pred):(n.pred+1+n.fixed.pred),]))
              X<-cbind(1, fixed.pred, pred)
              #V<-diag(rep(vy, times=N))+na.exclude(me.response)+diag(as.numeric(me.pred%*%(beta1[(2+n.fixed.pred):(n.pred+n.fixed.pred+1),]*beta1[(2+n.fixed.pred):(n.pred+n.fixed.pred+1),])))-diag(as.numeric(me.cov%*%(2*beta1[(2+n.fixed.pred):(n.pred+n.fixed.pred+1),]))) + diag(as.numeric(me.fixed.pred%*%(beta1[2:(length(beta1)-n.pred),]*beta1[2:(length(beta1)-n.pred),])))-diag(as.numeric(me.fixed.cov%*%(2*beta1[2:(length(beta1)-n.pred),])));

              V<-diag(rep(vy, times=N))+na.exclude(me.response)+ obs_var_con-diag(as.numeric(me.cov%*%(2*beta1[(2+n.fixed.pred):(n.pred+n.fixed.pred+1),]))) -diag(as.numeric(me.fixed.cov%*%(2*beta1[2:(length(beta1)-n.pred),])));

            }
            else
            {
              if(ultrametric==TRUE)
                s1<-as.numeric(s.X%*%(beta1[(2+n.fixed.pred):(n.pred+n.fixed.pred+1),]*beta1[(2+n.fixed.pred):(n.pred+n.fixed.pred+1),]))
              else
                s1<-as.numeric(s.X%*%(beta1[(4+n.fixed.pred):(n.pred+n.fixed.pred+3),]*beta1[(4+n.fixed.pred):(n.pred+n.fixed.pred+3),]))

              for(p in 1:N)
              {
                for(q in 1:N)
                {
                  if(ta[q,p]==0)num.prob[q,p]=1 else num.prob[q,p]=(1-exp(-a*ta[q,p]))/(a*ta[q,p])
                }
              }
              cm1<-(s1/(2*a)+vy)*(1-exp(-2*a*ta))*exp(-a*tij)
              for(p in 1:N)
              {
                for(q in 1:N)
                {
                  cm2[p,q]<-(((1-exp(-a*T[p]))/(a*T[p]))*((1-exp(-a*T[q]))/(a*T[q]))-(exp(-a*tia[p, q])*(1-exp(-a*T[p]))/(a*T[q])+ exp(-a*tja[p, q])*(1-exp(-a*T[p]))/(a*T[p]))*(num.prob[p,q]))
                }
              }
              if(ultrametric==TRUE)
              {
                X<-cbind(1, fixed.pred, (1-(1-exp(-a*T))/(a*T))*pred)
                mv<-diag(rowSums(matrix(data=as.numeric(me.pred)*t(kronecker(beta1[(2+n.fixed.pred):(n.pred+n.fixed.pred+1), ], (1-(1-exp(-a*T))/(a*T)))^2), ncol=n.pred)))
                mcov<-diag(rowSums(matrix(data=as.numeric(me.cov)*t(kronecker(2*beta1[(2+n.fixed.pred):(n.pred+n.fixed.pred+1),], (1-(1-exp(-a*T))/(a*T)))), ncol=n.pred)))

                #V<-cm1+(s1*ta*cm2)+me.response+mv-mcov + diag(as.numeric(me.fixed.pred%*%(beta1[2:(length(beta1)-n.pred),]*beta1[2:(length(beta1)-n.pred),])))-diag(as.numeric(me.fixed.cov%*%(2*beta1[2:(length(beta1)-n.pred),])));

                V<-cm1+(s1*ta*cm2)+na.exclude(me.response)+ obs_var_con-mcov -diag(as.numeric(me.fixed.cov%*%(2*beta1[2:(length(beta1)-n.pred),])));

              }
              else
              {
                nu.X<-cbind(1-exp(-a*T), 1-exp(-a*T)-(1-(1-exp(-a*T))/(a*T)), exp(-a*T), fixed.pred, (1-(1-exp(-a*T))/(a*T))*pred)
                mv<-diag(rowSums(matrix(data=as.numeric(me.pred)*t(kronecker(beta1[(4+n.fixed.pred):(n.pred+n.fixed.pred+3), ], (1-(1-exp(-a*T))/(a*T)))^2), ncol=n.pred)))

                mcov<-diag(rowSums(matrix(data=as.numeric(me.cov)*t(kronecker(2*beta1[(4+n.fixed.pred):(n.pred+n.fixed.pred+3),], (1-(1-exp(-a*T))/(a*T)))), ncol=n.pred)))

                #V<-cm1+(s1*ta*cm2)+me.response+mv-mcov + diag(as.numeric(me.fixed.pred%*%(beta1[4:(length(beta1)-n.pred),]*beta1[4:(length(beta1)-n.pred),])))-diag(as.numeric(me.fixed.cov%*%(2*beta1[4:(length(beta1)-n.pred),])));

                V<-cm1+(s1*ta*cm2)+na.exclude(me.response)+ obs_var_con- mcov -diag(as.numeric(me.fixed.cov%*%(2*beta1[4:(length(beta1)-n.pred),])));

              }
            } # END OF ELSE CONDITION FOR HALF-LIFE = 0

            # INTERMEDIATE ESTIMATION OF OPTIMAL REGRESSION #

            V.inverse<-solve(V)
            if(half_life_values[i]==0)
            {
              beta.i<-pseudoinverse(t(X)%*%V.inverse%*%X)%*%(t(X)%*%V.inverse%*%Y)
              test<-matrix(nrow=(n.pred+n.fixed.pred+1))
              for(f in 1:(n.pred+n.fixed.pred+1))
              {
                if(abs(as.numeric(beta.i[f]-beta1[f]))<=convergence) test[f]=0 else test[f]=1
              }
              if(sum(test)==0) break
              con.count=con.count+1
              if(con.count >= 50)
              {
                message("Warning, estimates did not converge after 50 iterations, last estimates printed out")
                break
              }

              beta1<-beta.i
            }
            else
            {
              if(ultrametric==TRUE)
              {
                beta.i<-pseudoinverse(t(X)%*%V.inverse%*%X)%*%(t(X)%*%V.inverse%*%Y)
                test<-matrix(nrow=(n.pred+n.fixed.pred+1))
                for(f in 1:(n.pred+n.fixed.pred+1))
                {
                  if(abs(as.numeric(beta.i[f]-beta1[f]))<=convergence) test[f]=0 else test[f]=1
                }
                if(sum(test)==0) break
                con.count=con.count+1
                if(con.count >= 50)
                {
                  message("Warning, estimates did not converge after 50 iterations, last estimates printed out")
                  break
                }

                beta1<-beta.i
              }
              else
              {
                beta.i<-pseudoinverse(t(nu.X)%*%V.inverse%*%nu.X)%*%(t(nu.X)%*%V.inverse%*%Y)
                test<-matrix(nrow=(n.pred+n.fixed.pred))
                for(f in 4:(n.pred+n.fixed.pred+3))
                {
                  if(abs(as.numeric(beta.i[f]-beta1[f]))<=convergence) test[(f-3)]=0 else test[(f-3)]=1
                }
                if(sum(test)==0) break
                con.count=con.count+1
                if(con.count >= 50)
                {
                  message("Warning, estimates did not converge after 50 iterations, last estimates printed out")
                  break
                }

                beta1<-beta.i
              }
            }                          # END OF HALF-LIFE = 0 CONDITION #
          }                            # END OF ITERATED GLS REPEAT LOOP #
          beta1<-beta.i

          ### END OF ITERATED GLS ESTIMATION FOR BETA #

          if(half_life_values[i]==0)
          {
            s1<-as.numeric(s.X%*%(beta1[(2+n.fixed.pred):(n.pred+n.fixed.pred+1),]*beta1[(2+n.fixed.pred):(n.pred+n.fixed.pred+1),]))
            X<-cbind(1, fixed.pred,pred)
            #V<-diag(rep(vy, times=N))+me.response+diag(as.numeric(me.pred%*%(beta1[(2+n.fixed.pred):(n.pred+n.fixed.pred++1),]*beta1[(2+n.fixed.pred):(n.pred+n.fixed.pred+1),])))-diag(as.numeric(me.cov%*%(2*beta1[(2+n.fixed.pred):(n.pred+n.fixed.pred+1),]))) + diag(as.numeric(me.fixed.pred%*%(beta1[2:(length(beta1)-n.pred),]*beta1[2:(length(beta1)-n.pred),])))-diag(as.numeric(me.fixed.cov%*%(2*beta1[2:(length(beta1)-n.pred),])))

            V<-diag(rep(vy, times=N))+na.exclude(me.response)+ obs_var_con-diag(as.numeric(me.cov%*%(2*beta1[(2+n.fixed.pred):(n.pred+n.fixed.pred+1),]))) -diag(as.numeric(me.fixed.cov%*%(2*beta1[2:(length(beta1)-n.pred),])))

            V.inverse<-solve(V)
            eY<-X%*%beta1
            resid<-Y-eY;
            gof[i, k] <- -N/2*log(2*pi)-0.5*log(det(V))-0.5*(t(resid) %*% V.inverse%*%resid);
          }
          else
          {
            if(ultrametric==TRUE)
              s1<-as.numeric(s.X%*%(beta1[(2+n.fixed.pred):(n.pred+1+n.fixed.pred),]*beta1[(2+n.fixed.pred):(n.pred+1+n.fixed.pred),]))
            else
              s1<-as.numeric(s.X%*%(beta1[(4+n.fixed.pred):(n.pred+3+n.fixed.pred),]*beta1[(4+n.fixed.pred):(n.pred+3+n.fixed.pred),]))
            for(p in 1:N)
            {
              for(q in 1:N)
              {
                if(ta[q,p]==0)num.prob[q,p]=1 else num.prob[q,p]=(1-exp(-a*ta[q,p]))/(a*ta[q,p]);
              }
            }
            cm1<-(s1/(2*a)+vy)*(1-exp(-2*a*ta))*exp(-a*tij);
            for(p in 1:N)
            {
              for(q in 1:N)
              {
                cm2[p,q]<-(((1-exp(-a*T[p]))/(a*T[p]))*((1-exp(-a*T[q]))/(a*T[q]))-(exp(-a*tia[p, q])*(1-exp(-a*T[p]))/(a*T[q])+ exp(-a*tja[p, q])*(1-exp(-a*T[p]))/(a*T[p]))*(num.prob[p,q]));
              }
            }
            if(ultrametric==TRUE)
            {
              X<-cbind(1, fixed.pred, (1-(1-exp(-a*T))/(a*T))*pred)
              mv<-diag(rowSums(matrix(data=as.numeric(me.pred)*t(kronecker(beta1[(2+n.fixed.pred):(n.pred+1+n.fixed.pred), ], (1-(1-exp(-a*T))/(a*T)))^2), ncol=n.pred)))
              mcov<-diag(rowSums(matrix(data=as.numeric(me.cov)*t(kronecker(2*beta1[(2+n.fixed.pred):(n.pred+1+n.fixed.pred),], (1-(1-exp(-a*T))/(a*T)))), ncol=n.pred)))
              V<-cm1+(s1*ta*cm2)+me.response+mv-mcov+ diag(as.numeric(me.fixed.pred%*%(beta1[2:(length(beta1)-n.pred),]*beta1[2:(length(beta1)-n.pred),])))-diag(as.numeric(me.fixed.cov%*%(2*beta1[2:(length(beta1)-n.pred),])));
            }
            else
            {
              nu.X<-cbind(1-exp(-a*T), 1-exp(-a*T)-(1-(1-exp(-a*T))/(a*T)), exp(-a*T), fixed.pred, (1-(1-exp(-a*T))/(a*T))*pred)
              mv<-diag(rowSums(matrix(data=as.numeric(me.pred)*t(kronecker(beta1[(4+n.fixed.pred):(n.pred+3+n.fixed.pred), ], (1-(1-exp(-a*T))/(a*T)))^2), ncol=n.pred)))
              mcov<-diag(rowSums(matrix(data=as.numeric(me.cov)*t(kronecker(2*beta1[(4+n.fixed.pred):(n.pred+3+n.fixed.pred),], (1-(1-exp(-a*T))/(a*T)))), ncol=n.pred)))

              V<-cm1+(s1*ta*cm2)+me.response+mv-mcov + diag(as.numeric(me.fixed.pred%*%(beta1[4:(length(beta1)-n.pred),]*beta1[4:(length(beta1)-n.pred),])))-diag(as.numeric(me.fixed.cov%*%(2*beta1[4:(length(beta1)-n.pred),])));

              V<-cm1+(s1*ta*cm2)+na.exclude(me.response)+ obs_var_con-mcov -diag(as.numeric(me.fixed.cov%*%(2*beta1[4:(length(beta1)-n.pred),])));


            }
            V.inverse<-solve(V)
            if(ultrametric==TRUE)
              eY<-X%*%beta1
            else
              eY<-nu.X%*%beta1
            resid<-Y-eY;
            gof[i, k] <- -N/2*log(2*pi)-0.5*log(det(V))-0.5*(t(resid) %*% V.inverse%*%resid);
          }  # END OF CONDITION FOR HALF-LIFE = 0 #

          print(as.numeric(round(cbind(if(a!=0)log(2)/a else 0.00, vy, gof[i,k], t(beta1)), 4)))

        }
      }



      x<-rev(half_life_values)
      y<-vy_values
      z<-gof;
      ml<-max(z);
      for(i in 1:length(half_life_values))
      {
        for(j in 1:length(vy_values))
        {
          if(gof[i,j]==ml){alpha.est=log(2)/half_life_values[i]; vy.est=vy_values[j]}
        }
      }
      for(i in 1:length(half_life_values))
      {
        for(j in 1:length(vy_values))
        {
          if(gof[i,j]<=ml-support)gof[i, j]=ml-support;
        }
      }
      gof=gof-ml



      # FINAL OPTIMAL REGRESSION USING BEST ALPHA AND VY ESTIMATES #

      if(alpha.est==Inf)
      {
        gls.beta1<-glsyx.beta1<- solve(t(x.ols)%*%x.ols)%*%(t(x.ols)%*%Y)
        con.count<-0 # counter to break loop in the event of non-convergence
        repeat
        {
          s1<-as.numeric(s.X%*%(gls.beta1[(2+n.fixed.pred):(n.pred+1+n.fixed.pred),]*gls.beta1[(2+n.fixed.pred):(n.pred+1+n.fixed.pred),]))
          X<-cbind(1, fixed.pred, pred)
          #V<-diag(rep(vy, times=N))+me.response+diag(as.numeric(me.pred%*%(gls.beta1[(2+n.fixed.pred):(n.pred+1+n.fixed.pred),]*gls.beta1[(2+n.fixed.pred):(n.pred+1+n.fixed.pred),])))-diag(as.numeric(me.cov%*%(2*gls.beta1[(2+n.fixed.pred):length(gls.beta1),]))) + diag(as.numeric(me.fixed.pred%*%(gls.beta1[2:(length(gls.beta1)-n.pred),]*gls.beta1[2:(length(gls.beta1)-n.pred),])))-diag(as.numeric(me.fixed.cov%*%(2*gls.beta1[2:(length(gls.beta1)-n.pred),])));

          V<-diag(rep(vy, times=N))+na.exclude(me.response)+ obs_var_con -diag(as.numeric(me.cov%*%(2*gls.beta1[(2+n.fixed.pred):length(gls.beta1),]))) -diag(as.numeric(me.fixed.cov%*%(2*gls.beta1[2:(length(gls.beta1)-n.pred),])));


          V.inverse<-solve(V)
          beta.i.var<-ev.beta.i.var<-pseudoinverse(t(X)%*%V.inverse%*%X)
          beta.i<-beta.i.var%*%(t(X)%*%V.inverse%*%Y)
          test<-matrix(nrow=(n.pred+n.fixed.pred+1))
          for(f in 1:(n.pred+1+n.fixed.pred))
          {
            if(abs(as.numeric(beta.i[f]-gls.beta1[f]))<=convergence) test[f]=0 else test[f]=1
          }
          if(sum(test)==0) break
          con.count=con.count+1
          if(con.count >= 50)
          {
            message("Warning, estimates did not converge after 50 iterations, last estimates printed out")
            break
          }
          gls.beta1<-glsyx.beta1<-beta.i
        }
        gls.beta1<-glsyx.beta1<-beta.i
        X<-cbind(1, fixed.pred,pred)
        #V<-diag(rep(vy, times=N))+me.response+diag(as.numeric(me.pred%*%(gls.beta1[(2+n.fixed.pred):(n.pred+1+n.fixed.pred),]*gls.beta1[(2+n.fixed.pred):(n.pred+1+n.fixed.pred),])))-diag(as.numeric(me.cov%*%(2*gls.beta1[(2+n.fixed.pred):length(gls.beta1),]))) + diag(as.numeric(me.fixed.pred%*%(gls.beta1[2:(length(gls.beta1)-n.pred),]*gls.beta1[2:(length(gls.beta1)-n.pred),])))-diag(as.numeric(me.fixed.cov%*%(2*gls.beta1[2:(length(gls.beta1)-n.pred),])));

        V<-diag(rep(vy, times=N))+na.exclude(me.response)+ obs_var_con -diag(as.numeric(me.cov%*%(2*gls.beta1[(2+n.fixed.pred):length(gls.beta1),]))) -diag(as.numeric(me.fixed.cov%*%(2*gls.beta1[2:(length(gls.beta1)-n.pred),])));


        pred.mean<-X%*%gls.beta1
        g.mean<-(t(rep(1, times=N))%*%solve(V)%*%Y)/sum(solve(V));
        sst<-t(Y-g.mean)%*% solve(V)%*%(Y-g.mean)
        sse<-t(Y-pred.mean)%*%solve(V)%*%(Y-pred.mean)
        r.squared<-(sst-sse)/sst


        ###### Start of Bias correction ######
        print(gls.beta1)


        adj<-matrix(data=0, ncol= 1+ n.pred+n.fixed.pred, nrow=N)


        #if(ultrametric==TRUE)


        for(i in 1:(n.fixed.pred))
        {
          adj[,(1+i)]<-mean(X[,(1+i)])
        }

        for(i in 1:(n.pred))
        {
          adj[,(1+n.fixed.pred +i)] <- as.numeric(sigma.X.estimate(pred[,i], me.pred[,i], topology, times)[1]);
        }


        if(ultrametric==TRUE) Vu<-Vu else Vu<-Vu[-(1:(2*N)),-(1:(2*N))]
        if(ultrametric==TRUE) Vd<-Vd else Vd<-Vd[-(1:(2*N)),-(1:(2*N))]

        correction<-matrix(Vu%*%pseudoinverse(Vd+Vu)%*%(c(X)-c(adj)),  ncol=ncol(X), nrow=nrow(X), byrow=F)



        bias_corr<-pseudoinverse(t(X)%*%V.inverse%*%X)%*%t(X)%*%V.inverse%*%correction

        m<-length(glsyx.beta1)



        corrected_betas<-solve(diag(1,m,m)-bias_corr)%*% glsyx.beta1


        ###### End of Bias correction ######


      }

      else
      {
        if(ultrametric==TRUE)
          gls.beta1<-solve(t(x.ols)%*%x.ols)%*%(t(x.ols)%*%Y)
        else
          gls.beta1<-rbind(0, 0, solve(t(x.ols)%*%x.ols)%*%(t(x.ols)%*%Y));
        con.count<-0;
        repeat
        {
          if(ultrametric==TRUE)
            s1<-as.numeric(s.X%*%(gls.beta1[(2+n.fixed.pred):(n.pred+1+n.fixed.pred),]*gls.beta1[(2+n.fixed.pred):(n.pred+1+n.fixed.pred),]))
          else
            s1<-as.numeric(s.X%*%(gls.beta1[(4+n.fixed.pred):(n.pred+3+n.fixed.pred),]*gls.beta1[(4+n.fixed.pred):(n.pred+3+n.fixed.pred),]))
          for(p in 1:N)
          {
            for(q in 1:N)
            {
              if(ta[q,p]==0)num.prob[q,p]=1 else num.prob[q,p]=(1-exp(-alpha.est*ta[q,p]))/(alpha.est*ta[q,p])
            }
          }
          cm1<-(s1/(2*alpha.est)+vy.est)*(1-exp(-2*alpha.est*ta))*exp(-alpha.est*tij)
          for(p in 1:N)
          {
            for(q in 1:N)
            {
              cm2[p,q]<-(((1-exp(-alpha.est*T[p]))/(alpha.est*T[p]))*((1-exp(-alpha.est*T[q]))/(alpha.est*T[q]))-(exp(-alpha.est*tia[p, q])*(1-exp(-alpha.est*T[p]))/(alpha.est*T[q])+ exp(-alpha.est*tja[p, q])*(1-exp(-alpha.est*T[p]))/(alpha.est*T[p]))*(num.prob[p,q]))
            }
          }
          if(ultrametric==TRUE)
          {
            X<-cbind(1, fixed.pred,(1-(1-exp(-alpha.est*T))/(alpha.est*T))*pred)
            mv<-diag(rowSums(matrix(data=as.numeric(me.pred)*t(kronecker(gls.beta1[(2+n.fixed.pred):(n.pred+1+n.fixed.pred), ], (1-(1-exp(-alpha.est*T))/(alpha.est*T)))^2), ncol=n.pred)))
            mcov<-diag(rowSums(matrix(data=as.numeric(me.cov)*t(kronecker(2*gls.beta1[(2+n.fixed.pred):length(gls.beta1),], (1-(1-exp(-alpha.est*T))/(alpha.est*T)))), ncol=n.pred)))

            #V<-cm1+(s1*ta*cm2)+me.response+mv-mcov+ diag(as.numeric(me.fixed.pred%*%(gls.beta1[2:(length(gls.beta1)-n.pred),]*gls.beta1[2:(length(gls.beta1)-n.pred),])))-diag(as.numeric(me.fixed.cov%*%(2*gls.beta1[2:(length(gls.beta1)-n.pred),])));
            V<-cm1+(s1*ta*cm2)+na.exclude(me.response)+ obs_var_con-mcov-diag(as.numeric(me.fixed.cov%*%(2*gls.beta1[2:(length(gls.beta1)-n.pred),])));
          }
          else
          {
            nu.X<-cbind(1-exp(-alpha.est*T), 1-exp(-alpha.est*T)-(1-(1-exp(-alpha.est*T))/(alpha.est*T)), exp(-alpha.est*T), fixed.pred, (1-(1-exp(-alpha.est*T))/(alpha.est*T))*pred)
            mv<-diag(rowSums(matrix(data=as.numeric(me.pred)*t(kronecker(gls.beta1[(4+n.fixed.pred):(n.pred+3+n.fixed.pred), ], (1-(1-exp(-alpha.est*T))/(alpha.est*T)))^2), ncol=n.pred)))

            mcov<-diag(rowSums(matrix(data=as.numeric(me.cov)*t(kronecker(2*gls.beta1[(4+n.fixed.pred):length(gls.beta1),], (1-(1-exp(-a*T))/(a*T)))), ncol=n.pred)))
            #V<-cm1+(s1*ta*cm2)+me.response+mv-mcov + diag(as.numeric(me.fixed.pred%*%(gls.beta1[4:(length(gls.beta1)-n.pred),]*gls.beta1[4:(length(gls.beta1)-n.pred),])))-diag(as.numeric(me.fixed.cov%*%(2*gls.beta1[4:(length(gls.beta1)-n.pred),])));
            V<-cm1+(s1*ta*cm2)+na.exclude(me.response)+ obs_var_con-mcov -diag(as.numeric(me.fixed.cov%*%(2*gls.beta1[4:(length(gls.beta1)-n.pred),])));


          }

          V.inverse<-solve(V)


          if(ultrametric==TRUE)
          {
            beta.i.var<-pseudoinverse(t(X)%*%V.inverse%*%X)
            beta.i<-beta.i.var%*%(t(X)%*%V.inverse%*%Y)
            test<-matrix(nrow=(n.pred+1+n.fixed.pred))
            for(f in 1:(n.pred+1+n.fixed.pred))
            {
              if(abs(as.numeric(beta.i[f]-gls.beta1[f]))<=convergence) test[f]=0 else test[f]=1
            }
            if(sum(test)==0) break
            con.count=con.count+1
            if(con.count >= 50)
            {
              message("Warning, estimates did not converge after 50 iterations, last estimates printed out")
              break
            }


            gls.beta1<-beta.i
          }
          else
          {
            beta.i.var<-pseudoinverse(t(nu.X)%*%V.inverse%*%nu.X)
            beta.i<-beta.i.var%*%(t(nu.X)%*%V.inverse%*%Y)
            test<-matrix(nrow=(n.pred))
            for(f in 4:(n.pred+3+n.fixed.pred))
            {
              if(abs(as.numeric(beta.i[f]-beta1[f]))<=convergence) test[(f-3)]=0 else test[(f-3)]=1
            }
            if(sum(test)==0) break
            con.count=con.count+1
            if(con.count >= 50)
            {
              message("Warning, estimates did not converge after 50 iterations, last estimates printed out")
              break
            }

            beta1<-beta.i
          }
        }

        # Update obs_var_con matrix after new beta estimates #
        mybiglist <- list()

        xx<-seq(from=1, to=length(Vu[,1]), by=N)
        mybiglist <- list()


        for (i in seq(from=1, to=nrow(gls.beta1), by=1)){
          for (j in seq(from=1, to=nrow(gls.beta1), by=1)) {
            tmp <- list(error_condition[xx[i]:(i*N),xx[j]:(j*N)]* gls.beta1[i]* gls.beta1[j])
            mybiglist[xx[i]+j] <- tmp
          }

        }


        mybiglist <-rmNullObs(mybiglist)
        obs_var_con<-Reduce('+', mybiglist)



        # END OF ITERATED GLS LOOP #


        # CODE FOR SST, SSE AND R-SQUARED #

        if(ultrametric==TRUE)
          gls.beta1<-beta.i
        else
        {
          gls.beta1<-beta.i
          ind.par<-matrix(data=0, nrow=N, ncol=4, dimnames=list(NULL, c("Bo", "Bi.Xia", "Yo", "Sum")))
          ind.par[,1]<-beta.i[1]*nu.X[,1]
          ind.par[,2]<-(beta.i[2]*nu.X[,2])
          ind.par[,3]<-beta.i[3]*nu.X[,3]
          ind.par[,4]<-ind.par[,1]+ind.par[,2]+ind.par[,3]
          mean.Bo=mean(ind.par[,4])
        }

        if(ultrametric==TRUE)
        {
          X<-cbind(1, fixed.pred,(1-(1-exp(-alpha.est*T))/(alpha.est*T))*pred)
          mv<-diag(rowSums(matrix(data=as.numeric(me.pred)*t(kronecker(gls.beta1[(2+n.fixed.pred):(n.pred+1+n.fixed.pred), ], (1-(1-exp(-alpha.est*T))/(alpha.est*T)))^2), ncol=n.pred)))
          mcov<-diag(rowSums(matrix(data=as.numeric(me.cov)*t(kronecker(2*gls.beta1[(2+n.fixed.pred):length(gls.beta1),], (1-(1-exp(-alpha.est*T))/(alpha.est*T)))), ncol=n.pred)))
          #V<-cm1+(s1*ta*cm2)+me.response+mv-mcov+ diag(as.numeric(me.fixed.pred%*%(beta1[2:(length(gls.beta1)-n.pred),]*gls.beta1[2:(length(gls.beta1)-n.pred),])))-diag(as.numeric(me.fixed.cov%*%(2*gls.beta1[2:(length(gls.beta1)-n.pred),])))
          V<-cm1+(s1*ta*cm2)+na.exclude(me.response)+ obs_var_con-mcov-diag(as.numeric(me.fixed.cov%*%(2*gls.beta1[2:(length(gls.beta1)-n.pred),])))


          pred.mean<-X%*%gls.beta1

        }
        else
        {
          nu.X<-cbind(1-exp(-alpha.est*T), 1-exp(-alpha.est*T)-(1-(1-exp(-alpha.est*T))/(alpha.est*T)), exp(-alpha.est*T),fixed.pred, (1-(1-exp(-alpha.est*T))/(alpha.est*T))*pred)
          mv<-diag(rowSums(matrix(data=as.numeric(me.pred)*t(kronecker(beta1[(4+n.fixed.pred):(n.pred+3+n.fixed.pred), ], (1-(1-exp(-alpha.est*T))/(alpha.est*T)))^2), ncol=n.pred)))
          mcov<-diag(rowSums(matrix(data=as.numeric(me.cov)*t(kronecker(2*beta1[(4+n.fixed.pred):length(beta1),], (1-(1-exp(-alpha.est*T))/(alpha.est*T)))), ncol=n.pred)))
          V<-cm1+(s1*ta*cm2)+na.exclude(me.response)+ obs_var_con-mcov-diag(as.numeric(me.fixed.cov%*%(2*gls.beta1[4:(length(gls.beta1)-n.pred),])))

          pred.mean<-nu.X%*%gls.beta1
        }

        g.mean<-(t(rep(1, times=N))%*%solve(V)%*%Y)/sum(solve(V));
        sst<-t(Y-g.mean)%*% solve(V)%*%(Y-g.mean)
        sse<-t(Y-pred.mean)%*%solve(V)%*%(Y-pred.mean)
        r.squared<-(sst-sse)/sst






        # FINAL EVOLUTIONARY REGRESSION USING BEST ALPHA AND VY ESTIMATES AND KNOWN VARIANCE MATRIX #


        if(ultrametric==TRUE)  s1<-as.numeric(s.X%*%(gls.beta1[(2+n.fixed.pred):(n.pred+1+n.fixed.pred),]*gls.beta1[(2+n.fixed.pred):(n.pred+1+n.fixed.pred),]))
        else s1<-as.numeric(s.X%*%(gls.beta1[(4+n.fixed.pred):(n.pred+3+n.fixed.pred),]*gls.beta1[(4+n.fixed.pred):(n.pred+3+n.fixed.pred),]));
        for(p in 1:N)
        {
          for(q in 1:N)
          {
            if(ta[q,p]==0)num.prob[q,p]=1 else num.prob[q,p]=(1-exp(-alpha.est*ta[q,p]))/(alpha.est*ta[q,p])
          }
        }
        cm1<-(s1/(2*alpha.est)+vy.est)*(1-exp(-2*alpha.est*ta))*exp(-alpha.est*tij)
        for(p in 1:N)
        {
          for(q in 1:N)
          {
            cm2[p,q]<-(((1-exp(-alpha.est*T[p]))/(alpha.est*T[p]))*((1-exp(-alpha.est*T[q]))/(alpha.est*T[q]))-(exp(-alpha.est*tia[p, q])*(1-exp(-alpha.est*T[p]))/(alpha.est*T[q])+ exp(-alpha.est*tja[p, q])*(1-exp(-alpha.est*T[p]))/(alpha.est*T[p]))*(num.prob[p,q]))
          }
        }

        if(ultrametric==TRUE)
          #V<-cm1+(s1*ta*cm2)+me.response+diag(as.numeric(me.pred%*%(gls.beta1[(2+n.fixed.pred):(n.pred+1+n.fixed.pred),]*gls.beta1[(2+n.fixed.pred):(n.pred+1+n.fixed.pred),])))-diag(as.numeric(me.cov%*%(2*gls.beta1[(2+n.fixed.pred):length(gls.beta1),]))) + diag(as.numeric(me.fixed.pred%*%(gls.beta1[2:(length(gls.beta1)-n.pred),]*gls.beta1[2:(length(gls.beta1)-n.pred),])))-diag(as.numeric(me.fixed.cov%*%(2*gls.beta1[2:(length(gls.beta1)-n.pred),])))

          V<-cm1+(s1*ta*cm2)+na.exclude(me.response)+ obs_var_con-diag(as.numeric(me.cov%*%(2*gls.beta1[(2+n.fixed.pred):length(gls.beta1),]))) -diag(as.numeric(me.fixed.cov%*%(2*gls.beta1[2:(length(gls.beta1)-n.pred),])))


        else
          #V<-cm1+(s1*ta*cm2)+me.response+diag(as.numeric(me.pred%*%(gls.beta1[(4+n.fixed.pred):(n.pred+3+n.fixed.pred),]*gls.beta1[(4+n.fixed.pred):(n.pred+3+n.fixed.pred),])))-diag(as.numeric(me.cov%*%(2*gls.beta1[(4+n.fixed.pred):length(gls.beta1),])))+ diag(as.numeric(me.fixed.pred%*%(gls.beta1[4:(length(gls.beta1)-n.pred),]*gls.beta1[4:(length(gls.beta1)-n.pred),])))-diag(as.numeric(me.fixed.cov%*%(2*gls.beta1[4:(length(gls.beta1)-n.pred),])))

          V<-cm1+(s1*ta*cm2)+na.exclude(me.response)+ obs_var_con-diag(as.numeric(me.cov%*%(2*gls.beta1[(4+n.fixed.pred):length(gls.beta1),])))-diag(as.numeric(me.fixed.cov%*%(2*gls.beta1[4:(length(gls.beta1)-n.pred),])))


        X1<-cbind(1, fixed.pred, pred)
        V.inverse<-solve(V)
        ev.beta.i.var<-pseudoinverse(t(X1)%*%V.inverse%*%X1)
        ev.beta.i<-ev.beta.i.var%*%(t(X1)%*%V.inverse%*%Y)
        glsyx.beta1<-ev.beta.i


        ###### Start of Bias correction ######



        adj<-matrix(data=0, ncol= 1+ n.pred+n.fixed.pred, nrow=N)


        #if(ultrametric==TRUE)


        for(i in 1:(n.fixed.pred))
        {
          adj[,(1+i)]<-mean(X[,(1+i)])
        }

        for(i in 1:(n.pred))
        {
          adj[,(1+n.fixed.pred +i)] <- as.numeric(sigma.X.estimate(pred[,i], me.pred[,i], topology, times)[1]);
        }


        if(ultrametric==TRUE) Vu<-Vu else Vu<-Vu[(length(Vu[1,])-(N*2)+1):length(Vu[1,]),(length(Vu[1,])-(N*2)+1):length(Vu[1,])]
        if(ultrametric==TRUE) Vd<-Vd else Vd<-Vd[(length(Vd[1,])-(N*2)+1):length(Vd[1,]),(length(Vd[1,])-(N*2)+1):length(Vd[1,])]



        correction<-matrix(Vu%*%pseudoinverse(Vd+Vu)%*%(c(X1)-c(adj)),  ncol=ncol(X1), nrow=nrow(X1), byrow=F)



        bias_corr<-pseudoinverse(t(X1)%*%V.inverse%*%X1)%*%t(X1)%*%V.inverse%*%correction

        m<-length(glsyx.beta1)



        corrected_betas<-solve(diag(1,m,m)-bias_corr)%*% glsyx.beta1


        ###### End of Bias correction ######

      }                                         # END OF HALFLIFE 0 CONDITION #



    } # END OF RANDOM AND FIXED COVARIATE REGRESSION ESTIMATION

  }# END OF FIXED AND RANDOM COVARIATE ANCOVA AND REGRESSION PARAMETER ESTIMATION
  #### END OF NEW CODE


  # PLOT THE SUPPORT SURFACE FOR HALF-LIVES AND VY


  if(length(half_life_values) > 1 && length(vy_values) > 1){
    z1<-gof
    for(i in 1:length(vy_values)){
      h.lives[,i]=rev(z1[,i])
    }
    z<-h.lives
    op <- par(bg = "white")

    plot.slouch.x <<- x
    plot.slouch.y <<- y
    plot.slouch.loglik <<- z

    persp(x, y, z, theta = plot.angle, phi = 30, expand = 0.5, col = "NA") ## plot.angle = 30 default
    persp(x, y, z, theta = plot.angle, phi = 30, expand = 0.5, col = "NA",
          ltheta = 120, shade = 0.75, ticktype = "detailed",
          xlab = "half-life", ylab = "vy", zlab = "log-likelihood") -> res
  }
  #plot.coord <- cbind(rep(x, length(y)), rep(rev(y), length(x)))




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


    optima[,1] = gls.beta1;
    optima[,2] = std;
    print(optima);message("");message("");
    print(corrected_beta_values);message("");

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


} # END OF MODEL FITTING FUNCTION



####### WEIGHT MATRIX FOR FACTORS, BORROWED AND MODIFIED FROM OUCH #########
weight.matrix<-function(alpha, topology, times, N, regime.specs, fixed.cov, intercept)
{
  if (alpha == Inf) alpha<-10000000000000000000
  N <- N
  reg <- set.of.regimes(topology, regime.specs)
  R <- length(reg)
  T <- times[terminal.twigs(topology)]
  ep <- epochs(topology, times, terminal.twigs(topology))
  beta <- regimes(topology, times, regime.specs, terminal.twigs(topology))
  W <- matrix(data = 0, nrow = N, ncol = R + 1, dimnames = list(c(),
                                                                c("Ya", as.character(set.of.regimes(topology, regime.specs)))))
  W[, 1] <- exp(-alpha * T)
  for (i in 1:N) {
    delta <- diff(exp(alpha * (ep[[i]] - T[i])))
    for (k in 1:R) {
      W[i, k + 1] <- -sum(delta * beta[[i + N * (k - 1)]])
    }
  }
  if (is.null(intercept))
    W <- W
  else {
    if (intercept == "root") {
      root.reg <- as.character(regime.specs[times == 0])
      nonroot.reg <- as.character(reg[reg != root.reg])
      int <- as.matrix(W[, 1] + W[, root.reg])
      colnames(int) = root.reg
      if (max(int[, 1]) <= 0.01)
        W <- W
      else {
        W2 <- cbind(int, W[, nonroot.reg])
        W <- W2
        if(max(abs(W[,1]-W[,2]))<= 0.01) W<-W[,(2:length(W[1,]))] #sjekk denne for treghet

        #W <- W2[, 2:length(W2[1,])]
      }
    }
    else W[, 1] <- intercept
  }
  if (max(W[, 1]) <= 0.01)
    W <- W[, -1];

  if(!is.null(fixed.cov))
  {
    fixed.pred<-data.frame(fixed.cov);
    n.fixed.pred<-length(fixed.pred[1,]);
    fixed.pred<-matrix(data=fixed.pred[!is.na(fixed.pred)], ncol=n.fixed.pred);
    W<-cbind(W, fixed.pred)
  }
  return(W)
}

np.regression<-function(response, me.response, predictor, me.predictor, convergence=NULL){
  if(is.null(convergence)) convergence=0.000001
  Y <- response[!is.na(response)];
  N <- length(Y);
  pred<-data.frame(predictor)
  me.pred<-data.frame(me.predictor)
  n.pred<-length(pred[1,])
  pred<-matrix(data=pred[!is.na(pred)], ncol=n.pred)
  me.pred<-matrix(data=me.pred[!is.na(me.pred)], ncol=n.pred)
  X<-cbind(1, pred)
  me1<-diag(me.response[!is.na(me.response)])
  x.ols<-cbind(1, pred)
  V1<-diag(rep(1, times=N))
  beta1<-solve(t(x.ols)%*%x.ols)%*%(t(x.ols)%*%Y)
  eY<-X%*%beta1;
  pred.mean<-X%*%beta1
  g.mean<-(t(rep(1, times=N))%*%solve(V1)%*%Y)/sum(solve(V1));
  sst<-t(Y-g.mean)%*% solve(V1)%*%(Y-g.mean)
  sse<-t(Y-pred.mean)%*%solve(V1)%*%(Y-pred.mean)
  sigma<-sse/(N-(n.pred+1))
  repeat{
    V<-diag(rep(sigma, times=N))+me1 + diag(as.numeric(me.pred%*%(beta1[2:(n.pred+1),]*beta1[2:(n.pred+1),])))
    V.inverse<-solve(V)
    beta.i<-pseudoinverse(t(X)%*%V.inverse%*%X)%*%(t(X)%*%V.inverse%*%Y)
    test<-matrix(nrow=(n.pred+1))
    for(f in 1:(n.pred+1))
    {
      if(abs(as.numeric(beta.i[f]-beta1[f]))<=convergence) test[f]=0 else test[f]=1
    }
    if(sum(test)==0) break
    beta1<-beta.i
  }
  beta1<-beta.i
  beta.i.var<-pseudoinverse(t(X)%*%V.inverse%*%X)
  eY<-X%*%beta1;
  pred.mean<-X%*%beta1
  g.mean<-(t(rep(1, times=N))%*%solve(V)%*%Y)/sum(solve(V));
  sst<-t(Y-g.mean)%*% solve(V)%*%(Y-g.mean)
  sse<-t(Y-pred.mean)%*%solve(V)%*%(Y-pred.mean)
  sigma<-sse/(N-(n.pred+1))
  r.squared<-(sst-sse)/sst

  sigma.sq.y<-(1/N)*(sum((Y-X%*%beta1)^2))
  log.like<-(-N/2)*log(2*pi)-(N/2)*log(sigma.sq.y)-sum((Y-X%*%beta1)^2/(2*sigma.sq.y))
  ml<-log.like
  reg.par<-matrix(data=0, nrow=(n.pred+1), ncol=2, dimnames=list(c("Intercept", if(n.pred==1) deparse(substitute(predictor)) else colnames(predictor)), c("Estimate", "Std. Error")))
  reg.par[,1]<-round(beta1,5);
  reg.par[,2]<-round(sqrt(diag(beta.i.var)),5)
  modfit<-matrix(data=0, nrow=7, ncol=1, dimnames=list(c("Support", "AIC", "AICc", "SIC", "r squared", "SST", "SSE"),("Value")))
  n.par=n.pred;
  modfit[1,1]=ml
  modfit[2,1]=-2*ml+2*(2+n.par)
  modfit[3,1]=modfit[2,1]+(2*(2+n.par)*((2+n.par)+1))/(N-(2+n.par)-1)
  modfit[4,1]=-2*ml+log(N)*(2+n.par)
  modfit[5,1]=r.squared*100
  modfit[6,1]=sst
  modfit[7,1]=sse
  message("REGRESSION PARAMETERS");message("");
  print(reg.par);message("");
  message("--------------------------------------------------");
  message("MODEL FIT");message("");
  print(modfit); message("");
  message("==================================================");
}

fitch<-function(tree.data, niche, deltran=FALSE, acctran=FALSE, root=NULL){
  niche<-as.factor(niche)
  anc.char<-as.character(tree.data$ancestor)
  anc.num<-as.numeric(anc.char)
  anc.num[1]<-0
  node.char<-as.character(tree.data$nodes)
  node.num<-as.numeric(node.char)
  niche.num<-as.numeric(niche)
  niches<-levels(niche)
  num.code<-1:length(levels(niche))
  translator<-cbind(niches, num.code)
  obs.char<-unique(niche.num[!is.na(niche.num)])
  N<-length(node.num)


  tree.matrix<-matrix(data=c(anc.num, node.num, niche.num), ncol=3,)
  colnames(tree.matrix)<-c("Ancestors", "Nodes", "States")
  Nint<-max(tree.matrix[,1])
  Nnodes<-max(tree.matrix[,2])



  downpass.node.states<-list()
  cost<-0

  # place tip character states in downpass list

  for(i in 1:Nnodes)
  {
    downpass.node.states[[i]]<-tree.matrix[i,3]
  }


  check.node.order<-tree.matrix[,2]-tree.matrix[,1]
  if(any(check.node.order<0))
    traversal<-2:Nint else traversal<-Nint:1

  # Fitch dwonpass algorithm for internal nodes

  for(i in traversal)
  {

    children<-which(tree.matrix[i,2]==tree.matrix[,1])


    if(length(children)==2)
    {
      child1.state<-downpass.node.states[[children[1]]];
      child2.state<-downpass.node.states[[children[2]]];
      downpass.node.states[[i]]<-intersect(child1.state, child2.state); 	 }

    # need to extend this for cases greater than tritomies;

    if(length(children) > 2)
    {
      child1.state<-downpass.node.states[[children[1]]];
      child2.state<-downpass.node.states[[children[2]]]; 	  child3.state<-downpass.node.states[[children[3]]];
      tmp<-intersect(child1.state, child2.state)
      downpass.node.states[[i]]<-intersect(tmp, child3.state); 	 }

    # if no intersect of states, use union and update costs for state change;

    if(length(downpass.node.states[[i]])==0)
    {
      if(length(children)==2)
      {
        downpass.node.states[[i]]<-union(child1.state, child2.state);
        cost = cost+1;
      }
      if(length(children)>2)
      {
        tmp<-union(child1.state, child2.state);
        downpass.node.states[[i]]<-union(tmp, child3.state);
        cost = cost+1;
      }
    }
  }

  # final traversal for root node if internal node order is reversed

  if(any(check.node.order<0))
  {
    children<-which(tree.matrix[1,2]==tree.matrix[,1])

    if(length(children)==2)
    {
      child1.state<-downpass.node.states[[children[1]]];
      child2.state<-downpass.node.states[[children[2]]];
      downpass.node.states[[1]]<-intersect(child1.state, child2.state); 	 }

    # need to extend this for any number of polytomies;

    if(length(children) > 2)
    {
      child1.state<-downpass.node.states[[children[1]]];
      child2.state<-downpass.node.states[[children[2]]]; 	  child3.state<-downpass.node.states[[children[3]]];
      downpass.node.states[[1]]<-intersect(c(child1.state, child2.state), child3.state); 	 }

    if(length(downpass.node.states[[1]])==0)
    {
      if(length(children)==2)
      {
        downpass.node.states[[1]]<-union(child1.state, child2.state);
        cost = cost+1;
      }
      if(length(children)>2)
      {
        tmp<-union(child1.state, child.state)
        downpass.node.states[[1]]<-union(tmp, child3.state);
        cost = cost+1;
      }
    }
  }  # end of reverse tree traversal with root done separately;


  # Fitch up pass algorithm
  # the downpass is not guaranteed to give the optimal
  # node states (as defined by minimum number of changes)
  # for this we need the slighly more complex Fitch up-pass
  # algorithm;

  if(any(check.node.order<0))
    traversal<-2:Nint else traversal<-Nint:2
  pre.traversal<-rev(traversal)


  # start with assuming that final state set for the root
  # is the same as the downpass state set for the root
  # note that if the root is ambiguous, we need to set it to a single state for the algorithm to work

  finalpass.node.states<-downpass.node.states



  if(!is.null(root)) finalpass.node.states[[1]]<-root

  if(length(finalpass.node.states[[1]]) >=2)
  {
    message("There is an ambiguity at the root, as given below")
    print(as.numeric(finalpass.node.states[[1]]))
    message("One needs to set this to one of the states using root = state in the function call before attempting the deltran or acctran reconstructions")
  }

  for(i in pre.traversal)
  {

    #parent<-which(tree.matrix[i,1]==tree.matrix[,2]);]
    parent<-tree.matrix[i,2]
    ancestor<-tree.matrix[i,1]
    children<-which(tree.matrix[,1]==tree.matrix[i,2]);

    # set final pass state for node i to the intersect of node i with its immediate ancestor's state

    finalpass.node.states[[parent]]<-intersect(downpass.node.states[[parent]], finalpass.node.states[[ancestor]]);

    if(setequal(finalpass.node.states[[parent]], finalpass.node.states[[ancestor]])==F)
    {
      # conditions for bifurcating

      if(length(children)==2)
      {
        child1.state<-downpass.node.states[[children[1]]];
        child2.state<-downpass.node.states[[children[2]]];
        if(length(intersect(child1.state, child2.state))!=0)
        {
          finalpass.node.states[[parent]]<-union(downpass.node.states[[parent]], (intersect(finalpass.node.states[[ancestor]], union(child1.state, child2.state))));
        }
        if(length(intersect(child1.state, child2.state))==0)
        {
          finalpass.node.states[[i]]<-union(downpass.node.states[[parent]], finalpass.node.states[[ancestor]])
        }

        if(deltran==TRUE){
          if(length(finalpass.node.states[[i]]) >=2)  finalpass.node.states[[i]]= finalpass.node.states[[ancestor]]}

        if(acctran==TRUE){
          if(length(finalpass.node.states[[i]]) >=2)  finalpass.node.states[[i]]= setdiff(finalpass.node.states[[parent]],finalpass.node.states[[ancestor]]) }


      }       # end of bifurcaring condition


      # condition for polytomy

      if(length(children)>2)
      {
        child1.state<-downpass.node.states[[children[1]]];
        child2.state<-downpass.node.states[[children[2]]];
        child3.state<-downpass.node.states[[children[3]]]

        if(length(intersect(child3.state, intersect(child1.state, child2.state)))!=0)
        {
          finalpass.node.states[[parent]]<-union(downpass.node.states[[parent]], (intersect(finalpass.node.states[[ancestor]], union(union(child3.state,  child1.state), child2.state))));
        }
        if(length(intersect(child3.state, intersect(child1.state, child2.state)))==0)
        {
          finalpass.node.states[[i]]<-union(downpass.node.states[[parent]], finalpass.node.states[[ancestor]])
        }

      }     # end of polytomy condition

    } # end of final parent state not equal ancestor state


  } # end of post tree traversal

  # enumerate all posible optimal state
  # reconstructions in dataframe (ouchtree) format
  # and reconstitue niches as characters;

  return(list(treematrix=tree.matrix, Final.states=finalpass.node.states, downpass.cost=cost, niche.code=translator))
}

make.states<-function(pars.object){
  x<-pars.object$Final.states
  n<-length(x)
  n.states<-length(pars.object$niche.code[,1])

  #which nodes are ambiguous

  count<-0

  ambig<-NA
  for(i in 1:n)
  {
    if(length(x[[i]])>=2)
    {
      count=count+1;
      ambig<-c(ambig,i)
    }
  }
  ambig=ambig[-1]
  n.ambig<-length(ambig)

  # choose first character state

  tmp<-matrix(data=NA, nrow=n, ncol=1)
  for(i in 1:n)
  {
    tmp[i,1]<-x[[i]][1]
    for(j in 1:n.states)
    {
      if (tmp[i,1]==j) tmp[i,1]<-pars.object$niche.code[,1][[j]]
    }
  }

  # encode ambiguous states as character ambiguous

  if(n.ambig!=0)
  {
    for(i in 1:n)
    {
      for(j in 1:n.ambig)
      {
        if(i==ambig[j]) tmp[i,1]="ambiguous"
      }
    }
  }

  return(as.factor(tmp))
}

n.ambig<-function(pars.states)
{
  states<-pars.states$Final.states
  N<-length(states)
  count<-0
  ambig<-NA
  for(i in 1:N)
  {
    if(length(states[[i]]) >=2)
    {
      count=count+1;
      ambig<-c(ambig, i)
    }

  }
  return(list(N.ambiguous=count, ambiguous.nodes=ambig[-1]))
}


`arrange.tree` <-
  function (root, topology) {
    k <- which(topology==root);
    n <- length(k);
    reltree <- rep(0,length(topology));
    reltree[root] <- 0.5;
    p <- list(NULL);
    if (n > 0) {
      m <- rep(0,n);
      for (j in 1:n) {
        p[[j]] <- arrange.tree(k[j],topology);
        m[j] <- length(which(p[[j]] != 0));
      }
      cm <- c(0,cumsum(m));
      for (j in 1:n) {
        reltree <- reltree + (cm[j]/sum(m))*(p[[j]] != 0) + (m[j]/sum(m))*p[[j]];
      }
    }
    return(reltree);
  }


`branch.times` <-
  function (topology, times) {
    term <- terminal.twigs(topology);
    N <- length(term);
    bt <- matrix(data=0,nrow=N,ncol=N);
    bt[1,1] <- times[term[1]];
    for (i in 2:N) {
      pedi <- pedigree(topology,term[i]);
      for (j in 1:(i-1)) {
        pedj <- pedigree(topology,term[j]);
        for (k in 1:length(pedi)) {
          if (any(pedj == pedi[k])) break;
        }
        bt[j,i] <- bt[i,j] <- times[pedi[k]];
      }
      bt[i,i] <- times[term[i]];
    }
    return(bt);
  }

`distance.matrix` <-
  function (topology, times) {
    term <- terminal.twigs(topology);
    N <- length(term);
    dm <- matrix(data=0,nrow=N,ncol=N);
    dm[1,1] <- 0;
    for (i in 2:N) {
      pedi <- pedigree(topology,term[i]);
      for (j in 1:(i-1)) {
        pedj <- pedigree(topology,term[j]);
        for (k in 1:length(pedi)) {
          if (any(pedj == pedi[k])) break;
        }
        dm[j,i] <- dm[i,j] <- (times[term[i]]-times[pedi[k]]) + (times[term[j]]-times[pedi[k]]);
      }
      dm[i,i] <- 0;
    }
    return(dm);
  }


`epochs` <-
  function (topology, times, term) {
    N <- length(term);
    e <- vector(length=N,mode="list");
    for (k in 1:N) {
      p <- pedigree(topology,term[k]);
      e[[k]] <- times[p];
    }
    return(e);
  }


`make.tree` <-
  function(n.tips, stretch=NULL){
    n.nodes<-n.tips*2-1
    int.nodes<-n.nodes-n.tips
    species<-c(rep(NA, times=int.nodes), letters[1:n.tips])
    generic<-c(0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13, 14, 14, 15, 15, 16, 16, 17, 17, 18, 18, 19, 19, 20, 20, 21, 21, 22, 22, 23, 23, 24, 24, 25, 25, 26, 26, 27, 27, 28, 28, 29, 29, 30, 30, 31, 31, 32, 32, 33, 33, 34, 34, 35, 35, 36, 36, 37, 37, 38, 38, 39, 39, 40, 40, 41, 41, 42, 42, 43, 43, 44, 44, 45, 45, 46, 46, 47, 47, 48, 48, 49, 49, 50, 50, 51, 51, 52, 52, 53, 53, 54, 54, 55, 55, 56, 56, 57, 57, 58, 58, 59, 59, 60, 60, 61, 61, 62, 62, 63, 63)
    if(n.tips==4)
    {
      ancestor<-generic[1:n.nodes]
      time<-c(0, 0.5, 0.5, 1, 1, 1, 1)
    }
    if(n.tips==6)
    {
      ancestor<-c(0, 1, 2, 2, 1, 3, 3, 4, 4, 5, 5)
      time<-c(0, 1/3, 2/3, 2/3, 2/3, rep(1, times=n.tips))
    }
    if(n.tips==8)
    {
      ancestor<-generic[1:n.nodes]
      time<-c(0, 1/3, 1/3, 2/3, 2/3, 2/3, 2/3, rep(1, times=n.tips))
    }
    if(n.tips==10)
    {
      ancestor<-c(0, 1, 1, 2, 3, 3, 4, 4, 2, generic[9:n.nodes])
      time<-c(0, 0.25, 0.5, 0.5, 0.75, 0.75, 0.75, 0.75, 0.75, rep(1, times=n.tips) )
    }
    if(n.tips==12)
    {
      ancestor<-c(0, 1, 2, 2, 1, generic[6:n.nodes])
      time<-c(0, 0.25, 0.5, 0.5, 0.5, 0.75, 0.75, 0.75, 0.75, 0.75, 0.75, rep(1, times=n.tips))
    }
    if(n.tips==14)
    {
      ancestor<-c(0, 1, 1, 2, 2, 3, 4, 4, 5, 5, 6, 6, 3, generic[14:n.nodes])
      time<-c(0, 0.25, 0.25, 0.5, 0.5, 0.5, 0.75, 0.75, 0.75, 0.75, 0.75, 0.75, 0.75, rep(1, times=n.tips))}
    if(n.tips==16)
    {
      ancestor<-generic[1:n.nodes]
      time<-c(0, 0.25, 0.25, rep(0.5, times=4), rep(0.75, times=8), rep(1, times=n.tips))
    }
    if(n.tips==18)
    {
      ancestor<-c(0, 1, 1, 3, 2, 2, 3, 4, 5, 5, 6, 6, 7, 7, 8, 8, 4, generic[18:n.nodes])
      time<-c(0, 0.2, 0.2, 0.4, rep(0.6, times=4), rep(0.8, times=9), rep(1, times=n.tips))
    }
    if(n.tips==20)
    {
      ancestor<-c(0, 1, 1, 2, 3, 3, 4, 4, 2, generic[10:n.nodes])
      time<-c(0, 0.2, 0.4, 0.4, rep(0.6, times=5), rep(0.8, times=10), rep(1, times=n.tips))
    }
    if(n.tips==32)
    {
      ancestor<-	generic[1:n.nodes]
      time<-c(0, 0.2, 0.2, rep(0.4, times=4), rep(0.6, times=8), rep(0.8, times=16), rep(1, times=n.tips))
      species<-c(rep(NA, times=int.nodes), 1:32)
    }
    if(n.tips==64)
    {
      ancestor<-generic[1:n.nodes]
      time<-c(0, 1/6, 1/6, rep(2/6, times=4), rep(0.5, times=8), rep(4/6, times=16), rep(5/6, times=32), rep(1, times=n.tips))
      species<-c(rep(NA, times=int.nodes), 1:64)
    }


    if(!is.null(stretch)) time<-time^stretch
    tree<-data.frame(ancestor, time, species)
    return(tree)
  }


`no.me.sigma.X.estimate` <-
  function (predictor, topology, times) {
    pt <- parse.tree(topology,times);
    n <- pt$N;
    v <- pt$bt;
    w <- matrix(data=1,nrow=pt$N,ncol=1);
    dat <- predictor[!is.na(predictor)];
    beta<-solve(t(w)%*%solve(v)%*%w)%*%(t(w)%*%solve(v)%*%dat)
    e<-dat-beta
    theta <- beta
    sigma <- sqrt((e %*% solve(v,e))/(n-1));
    dim(sigma) <- 1;
    return(list(as.numeric(theta), as.numeric(sigma)));
  }

`np.regression` <-
  function(response, me.response, predictor, me.predictor, convergence=NULL){
    if(is.null(convergence)) convergence=0.000001
    Y <- response[!is.na(response)];
    N <- length(Y);
    pred<-data.frame(predictor)
    me.pred<-data.frame(me.predictor)
    n.pred<-length(pred[1,])
    pred<-matrix(data=pred[!is.na(pred)], ncol=n.pred)
    me.pred<-matrix(data=me.pred[!is.na(me.pred)], ncol=n.pred)
    X<-cbind(1, pred)
    me1<-diag(me.response[!is.na(me.response)])
    x.ols<-cbind(1, pred)
    V1<-diag(rep(1, times=N))
    beta1<-solve(t(x.ols)%*%x.ols)%*%(t(x.ols)%*%Y)
    eY<-X%*%beta1;
    pred.mean<-X%*%beta1
    g.mean<-(t(rep(1, times=N))%*%solve(V1)%*%Y)/sum(solve(V1));
    sst<-t(Y-g.mean)%*% solve(V1)%*%(Y-g.mean)
    sse<-t(Y-pred.mean)%*%solve(V1)%*%(Y-pred.mean)
    sigma<-sse/(N-(n.pred+1))
    repeat{
      V<-diag(rep(sigma, times=N))+me1 + diag(as.numeric(me.pred%*%(beta1[2:(n.pred+1),]*beta1[2:(n.pred+1),])))
      V.inverse<-solve(V)
      beta.i<-pseudoinverse(t(X)%*%V.inverse%*%X)%*%(t(X)%*%V.inverse%*%Y)
      test<-matrix(nrow=(n.pred+1))
      for(f in 1:(n.pred+1))
      {
        if(abs(as.numeric(beta.i[f]-beta1[f]))<=convergence) test[f]=0 else test[f]=1
      }
      if(sum(test)==0) break
      beta1<-beta.i
    }
    beta1<-beta.i
    beta.i.var<-pseudoinverse(t(X)%*%V.inverse%*%X)
    eY<-X%*%beta1;
    pred.mean<-X%*%beta1
    g.mean<-(t(rep(1, times=N))%*%solve(V)%*%Y)/sum(solve(V));
    sst<-t(Y-g.mean)%*% solve(V)%*%(Y-g.mean)
    sse<-t(Y-pred.mean)%*%solve(V)%*%(Y-pred.mean)
    sigma<-sse/(N-(n.pred+1))
    r.squared<-(sst-sse)/sst

    sigma.sq.y<-(1/N)*(sum((Y-X%*%beta1)^2))
    log.like<-(-N/2)*log(2*pi)-(N/2)*log(sigma.sq.y)-sum((Y-X%*%beta1)^2/(2*sigma.sq.y))
    ml<-log.like
    reg.par<-matrix(data=0, nrow=(n.pred+1), ncol=2, dimnames=list(c("Intercept", if(n.pred==1) deparse(substitute(predictor)) else colnames(predictor)), c("Estimate", "Std. Error")))
    reg.par[,1]<-round(beta1,5);
    reg.par[,2]<-round(sqrt(diag(beta.i.var)),5)
    modfit<-matrix(data=0, nrow=7, ncol=1, dimnames=list(c("Support", "AIC", "AICc", "SIC", "r squared", "SST", "SSE"),("Value")))
    n.par=n.pred;
    modfit[1,1]=ml
    modfit[2,1]=-2*ml+2*(2+n.par)
    modfit[3,1]=modfit[2,1]+(2*(2+n.par)*((2+n.par)+1))/(N-(2+n.par)-1)
    modfit[4,1]=-2*ml+log(N)*(2+n.par)
    modfit[5,1]=r.squared*100
    modfit[6,1]=sst
    modfit[7,1]=sse
    message("REGRESSION PARAMETERS");message("");
    print(reg.par);message("");
    message("--------------------------------------------------");
    message("MODEL FIT");message("");
    print(modfit); message("");
    message("==================================================");
  }

`oubm.sim` <-
  function(Xo, Yo, halflife, vx, vy, increment, b0, b1, topology, times, specs, n){
    alpha<-log(2)/halflife
    s.X<-sqrt(vx)
    s.y<-sqrt((vy*2*alpha))
    phylogeny<-cbind(topology, times, specs)
    N<-max(topology)  # N internal nodes
    int.node<-rep(NA, times=N)
    p<-matrix(data=0, nrow=length(specs), ncol=n)
    r<-matrix(data=0, nrow=length(specs), ncol=n)
    for(i in 1:n){
      z<-phylogeny.evolution(Xo, Yo, alpha, s.X, s.y, increment, b0, b1, topology, times)
      x.values<-z[,1]
      y.values<-z[,2]
      p[,i]<-as.numeric(c(int.node, x.values))
      r[,i]<-as.numeric(c(int.node, y.values))
    }
    sim.dat<-data.frame(cbind(phylogeny, p, r))
    x.head<-paste(rep("predictor_", times=n), 1:n, sep="")
    y.head<-paste(rep("response_", times=n), 1:n, sep="")
    names(sim.dat)<-c("ancestor", "time", "species", x.head, y.head)
    return(sim.dat)
  }

`ouch2slouch` <-
  function(tree)
  {

    if (!inherits(tree, "ouchtree"))
      stop(sQuote("tree"), " must be of class ", sQuote("ouchtree"))

    N<-length(tree@nodes)
    tmp<-as(tree, "data.frame")
    tmp$ancestors<-as.character(tmp$ancestors)
    tmp$ancestors<-as.numeric(tmp$ancestors)
    tmp$times<-as.character(tmp$times)
    tmp$times<-as.numeric(tmp$times)
    tmp$nodes<-as.character(tmp$nodes)
    tmp$nodes<-as.numeric(tmp$nodes)
    tmp$ancestors[1]<-0;
    slouch_node<-1:N
    ancestor<-rep(NA, times=N)

    for(i in 1:N)
    {
      if(tmp$labels[i]=="") tmp$labels[i]= NA;
    }

    rownames(tmp)<-1:nrow(tmp)
    names(tmp)=c("nodes", "ancestor", "time", "species")
    return(tmp)
  }

`parse.tree` <-
  function(topology, times){
    term<-terminal.twigs(topology);
    N <- length(term);
    dm<-distance.matrix(topology, times);
    bt<-branch.times(topology, times);
    pt<-list(N=N, term=term, dm=dm, bt=bt);
    return(pt);
  }

`pedigree` <-
  function (topology, k) {
    p <- k;
    k <- topology[k];
    while (k != 0) {
      p <- c(p, k);
      k <- topology[k];
    }
    return(p);
  }


`phylogeny.evolution` <-
  function(Xo, Yo, alpha, s.X, s.y, increment, b0, b1, topology, times){
    branch<-matrix(data=0, nrow=length(times), ncol=2);
    branch[1, 1]<-Xo;
    branch[1, 2]<-Yo;
    sx<-sqrt(increment)*s.X
    sy<-sqrt(increment)*s.y
    for(i in 2:length(times)){
      t<-times[i]-times[topology[i]];
      m<-round(t/increment)                   # note that some discretization error is introduced here, but it shouldnt matter as long as increment is small enough;
      X<-branch[topology[i], 1];
      Y<-branch[topology[i], 2];
      for (j in 1:m){
        Y<-Y+(-alpha*(Y-(b0+b1*(X)))*increment)+rnorm(1, mean=0, sd=sy);
        X<-X+rnorm(1, mean=0, sd=sx);
      }
      branch[i,1]<-X
      branch[i,2]<-Y
    }
    tip.specs=length(topology)-max(topology);
    branch1<-branch[tip.specs:length(topology),]
    return(branch1)
  }


pseudoinverse <-
  function (m){
    msvd <- svd(m)
    if (length(msvd$d) == 0) {
      return(array(0, dim(m)[2:1]))
    }
    else {
      return(msvd$v %*% (1/msvd$d * t(msvd$u)))
    }
  }

`regimes` <-
  function (topology, times, regime.specs, term) {
    N <- length(term);
    reg <- set.of.regimes(topology,regime.specs);
    R <- length(reg);
    beta <- vector(R*N, mode="list");
    for (i in 1:N) {
      for (k in 1:R) {
        p <- pedigree(topology, term[i]);
        n <- length(p);
        beta[[i + N*(k-1)]] <- as.integer(regime.specs[p[1:(n-1)]] == reg[k]);
      }
    }
    return(beta);
  }


`set.of.regimes` <-
  function (topology, regime.specs) {
    n <- length(regime.specs);
    id <- seq(1,n)[topology > 0];       # find all non-root nodes
    reg <- sort(unique(regime.specs[id]));
    return(reg);
  }

`sigma.X.estimate` <-
  function (predictor,me.predictor, topology, times) {
    pt <- parse.tree(topology,times);
    n <- pt$N;
    v <- pt$bt;
    w <- matrix(data=1,nrow=pt$N,ncol=1);
    me<-diag(me.predictor[!is.na(me.predictor)])
    dat <- predictor[!is.na(predictor)];
    beta<-solve(t(w)%*%solve(v)%*%w)%*%(t(w)%*%solve(v)%*%dat)
    e<-dat-beta
    sigma<-as.numeric((t(e)%*%solve(v)%*%e)/(n-1))
    repeat{
      beta<-solve(t(w)%*%solve(v + me/sigma)%*%w)%*%(t(w)%*%solve(v + me/sigma)%*%dat)
      e<-dat-beta
      theta <- beta
      sigma1<-(t(e)%*%solve(v +me/sigma)%*%e)/(n-1)
      if(abs(as.numeric(sigma1)-sigma) <= 0.0000001*sigma) break
      sigma <- as.numeric(sigma1)
    }
    return(list(as.numeric(theta), as.numeric(sigma)));
  }

slouchtree.plot <-function (topology, times, names = NULL, regimes = NULL, cex = NULL, lwd=NULL, reg.col=NULL) {
  if(is.null(cex)) cex<-1;
  if(is.null(lwd)) lwd<-1;
  rx <- range(times);
  rxd <- 0.1*diff(rx);

  if (is.null(regimes))
    regimes <- factor(rep(1,length(topology)));

  levs <- levels(as.factor(regimes));
  palette <- rainbow(length(levs));

  for (r in 1:length(levs)) {
    y <- tree.layout(topology);
    x <- times;
    f <- which(topology > 0 & regimes == levs[r]);
    pp <- topology[f];
    X <- array(data=c(x[f], x[pp], rep(NA,length(f))),dim=c(length(f),3));
    Y <- array(data=c(y[f], y[pp], rep(NA,length(f))),dim=c(length(f),3));
    oz <- array(data=1,dim=c(2,1));
    X <- kronecker(t(X),oz);
    Y <- kronecker(t(Y),oz);
    X <- X[2:length(X)];
    Y <- Y[1:(length(Y)-1)];
    if(!is.null(regimes))
    {if(is.null(reg.col))
      C <- rep(palette[r],length(X))
    }
    {if(!is.null(reg.col))
      C <- rep(reg.col[r],length(X))
    }
    if (r > 1) par(new=TRUE);
    par(yaxt='n')
    par(bty="n")
    par(font="2")

    plot(X,Y,type='l',col=C,lwd=lwd,xlab='time',ylab='',xlim = rx + c(-rxd,rxd),ylim=c(0,1));
    if (!is.null(names))
      text(X[seq(1,length(X),6)],Y[seq(1,length(Y),6)],names[f],pos=4, cex=cex);
  }
  par(yaxt="s") #reset graphic parameter to default
  par(bty="o")
  par(font="1")
}


`terminal.twigs` <-
  function(topology){
    n<-length(topology);
    return(seq(max(topology)+1, n));
  }

`tree.layout` <-
  function (topology) {
    root <- which(topology==0);
    return(arrange.tree(root,topology));
  }

`tsia` <-
  function (topology, times) {
    term <- terminal.twigs(topology);
    N <- length(term);
    t.ia <- matrix(data=0,nrow=N,ncol=N);
    t.ia[1,1]=0
    for (i in 2:N) {
      pedi <- pedigree(topology,term[i]);
      for (j in 1:(i-1)) {
        pedj <- pedigree(topology,term[j]);
        for (k in 1:length(pedi)) {
          if (any(pedj == pedi[k])) break;
        }
        t.ia[j,i] <- t.ia[i,j] <- (times[term[i]]-times[pedi[k]]);
      }
      t.ia[i,i] <- 0;
    }
    return(t.ia);
  }

`tsja` <-
  function (topology, times) {
    term <- terminal.twigs(topology);
    N <- length(term);
    t.ja <- matrix(data=0,nrow=N,ncol=N);
    t.ja[1,1]=0
    for (i in 2:N) {
      pedi <- pedigree(topology,term[i]);
      for (j in 1:(i-1)) {
        pedj <- pedigree(topology,term[j]);
        for (k in 1:length(pedi)) {
          if (any(pedj == pedi[k])) break;
        }
        t.ja[j,i] <- t.ja[i,j] <- (times[term[j]]-times[pedi[k]]);
      }
      t.ja[i,i] <- 0;
    }
    return(t.ja);
  }


make.slouch.data<-function(tree, data)
{
  if (inherits(tree, "phylo"))
  {
    tree<-ape2ouch(tree)
    tree<-ouch2slouch(tree)
  }
  if (inherits(tree, "ouchtree"))
  {
    tree<-ouch2slouch(tree)
  }

  N.col.tree<-length(tree[1,]);
  N.row.tree<-length(tree[,1]);
  N.col.data<-length(data[1,]);
  N.int<-length(tree$species[is.na(tree$species)]);
  int.tree<-tree[1:N.int,];
  tip.tree<-tree[(N.int+1):N.row.tree, ];
  tip.dat<-merge(tip.tree, data, by.x="species", by.y="species", sort=FALSE);
  x<-matrix(NA, ncol=(N.col.data-1), nrow=N.int);
  int<-data.frame(int.tree$species, int.tree$nodes, int.tree$ancestor, int.tree$time, x);
  names(int)<-names(tip.dat);
  slouch.dat<-rbind(int, tip.dat);

  return(slouch.dat)
}

is.NullOb <- function(x) is.null(x) | all(sapply(x, is.null))


rmNullObs <- function(x)
{
  x <- Filter(Negate(is.NullOb), x)
  lapply(x, function(x) if (is.list(x)) rmNullObs(x) else x)
}

###############################################################

sapply_pb <- function(X, FUN, ...)
{
  env <- environment()
  pb_Total <- length(X)
  counter <- 0
  pb <- txtProgressBar(min = 0, max = pb_Total, style = 3)

  wrapper <- function(...){
    curVal <- get("counter", envir = env)
    assign("counter", curVal +1 ,envir=env)
    setTxtProgressBar(get("pb", envir=env), curVal +1)
    FUN(...)
  }
  res <- sapply(X, wrapper, ...)
  close(pb)
  res
}

