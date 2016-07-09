sup.rReg <- function(hl_vy, N, me.response, ta, tij, T, topology, times, regime.specs, model.type, ultrametric, Y, fixed.cov, pred, xx, beta1, error_condition, s.X, n.pred, num.prob, tia, tja, cm2, me.pred, me.cov, convergence, n.fixed) {


  hl <- hl_vy[1]; vy <- hl_vy[2]
  if(hl==0)
  {
    x.ols<-cbind(1, pred)
    beta1<-solve(t(x.ols)%*%x.ols)%*%(t(x.ols)%*%Y)
  }
  else
  {
    a <- log(2)/hl
    x.ols<-cbind(1, pred)
    if(ultrametric==TRUE)
      beta1<-solve(t(x.ols)%*%x.ols)%*%(t(x.ols)%*%Y)
    else
      beta1<-rbind(0, 0, solve(t(x.ols)%*%x.ols)%*%(t(x.ols)%*%Y))

    obs_var_con <-matrix(0, nrow=N, ncol=N)

    for (e in seq(from=1, to=ncol(x.ols), by=1)){
      for (j in seq(from=1, to=ncol(x.ols), by=1)) {
        tmp<-error_condition[xx[e]:(e*N),xx[j]:(j*N)]*beta1[e]*beta1[j]
        obs_var_con <-obs_var_con + tmp
      }

    }

  }

  ### CODE FOR ESTIMATING BETA USING ITERATED GLS ###
  con.count<-0;  # Counter for loop break if Beta's dont converge #
  repeat
  {
    if(hl==0)
    {
      a<-Inf
      s1<-as.numeric(s.X%*%(beta1[2:(n.pred+1),]*beta1[2:(n.pred+1),]))
      X<-cbind(1, pred)

      V<-diag(rep(vy, times=N))+ na.exclude(me.response) + obs_var_con-diag(as.vector(na.exclude(me.cov%*%(2*beta1[(n.fixed+1):length(beta1),]))));


    }
    else
    {
      if(ultrametric==TRUE)
        s1<-as.numeric(s.X%*%(beta1[2:(n.pred+1),]*beta1[2:(n.pred+1),]))
      else
        s1<-as.numeric(s.X%*%(beta1[4:(n.pred+3),]*beta1[4:(n.pred+3),]))

      obs_var_con <-matrix(0, nrow=N, ncol=N)

      for (e in seq(from=1, to=ncol(x.ols), by=1)){
        for (j in seq(from=1, to=ncol(x.ols), by=1)) {
          tmp<-error_condition[xx[e]:(e*N),xx[j]:(j*N)]*(beta1[e]*(1-(1-exp(-a*T))/(a*T)))*(beta1[j]*(1-(1-exp(-a*T))/(a*T)))
          obs_var_con <-obs_var_con + tmp
        }

      }

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
        X<-cbind(1, (1-(1-exp(-a*T))/(a*T))*pred)
        mv<-diag(rowSums(matrix(data=as.numeric(me.pred)*t(kronecker(beta1[2:(n.pred+1), ], (1-(1-exp(-a*T))/(a*T)))^2), ncol=n.pred)))
        mcov<-diag(rowSums(matrix(data=as.numeric(me.cov)*t(kronecker(2*beta1[2:(n.pred+1),], (1-(1-exp(-a*T))/(a*T)))), ncol=n.pred)))

        V<-cm1+(s1*ta*cm2)+na.exclude(me.response)+ obs_var_con -mcov

      }
      else
      {
        nu.X<-cbind(1-exp(-a*T), 1-exp(-a*T)-(1-(1-exp(-a*T))/(a*T)), exp(-a*T), (1-(1-exp(-a*T))/(a*T))*pred)
        mv<-diag(rowSums(matrix(data=as.numeric(me.pred)*t(kronecker(beta1[4:(n.pred+3), ], (1-(1-exp(-a*T))/(a*T)))^2), ncol=n.pred)))

        mcov<-diag(rowSums(matrix(data=as.numeric(me.cov)*t(kronecker(2*beta1[4:(n.pred+3),], (1-(1-exp(-a*T))/(a*T)))), ncol=n.pred)))

        #V<-cm1+(s1*ta*cm2)+na.exclude(me.response)+mv-mcov
        V<-cm1+(s1*ta*cm2)+na.exclude(me.response)+ obs_var_con-mcov


      }
    } # END OF ELSE CONDITION FOR HALF-LIFE = 0

    # INTERMEDIATE ESTIMATION OF OPTIMAL REGRESSION #
    V.inverse<-solve(V)

    if(hl==0)
    {
      beta.i<-pseudoinverse(t(X)%*%V.inverse%*%X)%*%(t(X)%*%V.inverse%*%Y)
      test<-matrix(nrow=(n.pred+1))
      for(f in 1:(n.pred+1))
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
        # print(beta.i)
        test<-matrix(nrow=(n.pred+1))
        for(f in 1:(n.pred+1))
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

        beta1<-beta.i;

      }
      else
      {
        beta.i<-pseudoinverse(t(nu.X)%*%V.inverse%*%nu.X)%*%(t(nu.X)%*%V.inverse%*%Y)
        ### PROBLEM: beta.i blir en vektor med NaN. Problemet oppstår ved Ultrametric = False når treet faktisk er ultrametrisk.
        test<-matrix(nrow=(n.pred))
        for(f in 4:(n.pred+3))
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


  ### END OF ITERATED GLS ESTIMATION FOR BETA #  #### NEw obs_var_con?

  if(hl==0)
  {
    s1<-as.numeric(s.X%*%(beta1[2:(n.pred+1),]*beta1[2:(n.pred+1),]))
    X<-cbind(1, pred)

    #V<-diag(rep(vy, times=N))+na.exclude(me.response)+diag(as.numeric(me.pred%*%(beta1[2:(n.pred+1),]*beta1[2:(n.pred+1),])))-diag(as.numeric(me.cov%*%(2*beta1[2:(n.pred+1),])))
    #V<-diag(rep(vy, times=N))+na.exclude(me.response)+((Vu-(Vu%*%pseudoinverse(Vu+Vd)%*%Vu))*(beta1[2:(n.pred+1),]*beta1[2:(n.pred+1),]))-diag(as.numeric(me.cov%*%(2*beta1[2:(n.pred+1),])))
    V<-diag(rep(vy, times=N))+ na.exclude(me.response) + obs_var_con-diag(as.vector(na.exclude(me.cov%*%(2*beta1[(n.fixed+1):length(beta1),]))));

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

    #gof[i, k] <- -N/2*log(2*pi)-0.5*log.det.V-0.5*(t(resid) %*% V.inverse%*%resid);
    sup1 <- -N/2*log(2*pi)-0.5*log.det.V-0.5*(t(resid) %*% V.inverse%*%resid)
  }
  else
  {

    if(ultrametric==TRUE)
      s1<-as.numeric(s.X%*%(beta1[2:(n.pred+1),]*beta1[2:(n.pred+1),]))
    else
      s1<-as.numeric(s.X%*%(beta1[4:(n.pred+3),]*beta1[4:(n.pred+3),]))


    obs_var_con <-matrix(0, nrow=N, ncol=N)

    for (e in seq(from=1, to=ncol(x.ols), by=1)){
      for (j in seq(from=1, to=ncol(x.ols), by=1)) {
        tmp<-error_condition[xx[e]:(e*N),xx[j]:(j*N)]*(beta1[e]*(1-(1-exp(-a*T))/(a*T)))*(beta1[j]*(1-(1-exp(-a*T))/(a*T)))
        obs_var_con <-obs_var_con + tmp
      }

    }

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


      X<-cbind(1, (1-(1-exp(-a*T))/(a*T))*pred)
      mv<-diag(rowSums(matrix(data=as.numeric(me.pred)*t(kronecker(beta1[2:(n.pred+1), ], (1-(1-exp(-a*T))/(a*T)))^2), ncol=n.pred)))
      mcov<-diag(rowSums(matrix(data=as.numeric(me.cov)*t(kronecker(2*beta1[2:(n.pred+1),], (1-(1-exp(-a*T))/(a*T)))), ncol=n.pred)))

      #V<-cm1+(s1*ta*cm2)+me.response+mv-mcov;
      V<-cm1+(s1*ta*cm2)+na.exclude(me.response)+ obs_var_con-mcov

    }
    else
    {
      nu.X<-cbind(1-exp(-a*T), 1-exp(-a*T)-(1-(1-exp(-a*T))/(a*T)), exp(-a*T), (1-(1-exp(-a*T))/(a*T))*pred)
      mv<-diag(rowSums(matrix(data=as.numeric(me.pred)*t(kronecker(beta1[4:(n.pred+3), ], (1-(1-exp(-a*T))/(a*T)))^2), ncol=n.pred)))
      mcov<-diag(rowSums(matrix(data=as.numeric(me.cov)*t(kronecker(2*beta1[4:(n.pred+3),], (1-(1-exp(-a*T))/(a*T)))), ncol=n.pred)))

      #V<-cm1+(s1*ta*cm2)+me.response+mv-mcov
      V<-cm1+(s1*ta*cm2)+na.exclude(me.response)+ obs_var_con-mcov

    }
    V.inverse<-solve(V)
    if(ultrametric==TRUE)
      eY<-X%*%beta1
    else
      eY<-nu.X%*%beta1
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
    # gof[i, k] <- -N/2*log(2*pi)-0.5*log.det.V-0.5*(t(resid) %*% V.inverse%*%resid);
    sup1 <- -N/2*log(2*pi)-0.5*log.det.V-0.5*(t(resid) %*% V.inverse%*%resid)

  }  # END OF CONDITION FOR HALF-LIFE = 0 #
  return(sup1)
  print(as.numeric(round(cbind(if(a!=0)log(2)/a else 0.00, vy, sup1, t(beta1)), 4))) # Will increasing the number of digits avoid problems when plotting grid?
}


