## Function to regress the mfReg

sup.mfReg <- function(hl_vy, N, me.response, ta, tij, T, topology, times, model.type, ultrametric, Y, fixed.cov, pred, xx, beta1, error_condition, s.X, n.pred, num.prob, tia, tja, cm2, me.pred, me.cov, convergence, n.fixed, fixed.pred, n.fixed.pred, obs_var_con, me.fixed.cov, me.fixed.pred){
  ## Split hl_vy vector
  hl <- hl_vy[1]; vy <- hl_vy[2]
  x.ols<-cbind(1,fixed.pred, pred)
  
  ## Initial OLS estimates
  if(hl==0)
  {
    a<-Inf
    beta1<-solve(t(x.ols)%*%x.ols)%*%(t(x.ols)%*%Y)
  }
  else
  {
    a <- log(2)/hl
    T.term <- T
    cm2 <- make.cm2(a,tia,tja,ta,N,T.term)
    if(ultrametric==TRUE)
      beta1<-solve(t(x.ols)%*%x.ols)%*%(t(x.ols)%*%Y)
    else
      beta1<-rbind(0, 0, solve(t(x.ols)%*%x.ols)%*%(t(x.ols)%*%Y))
  }
  
  ### CODE FOR ESTIMATING BETA USING ITERATED GLS ###
  con.count<-0;  # Counter for loop break if Beta's dont converge #
  repeat
  {
    if(hl==0)
    {
      s1<-as.numeric(s.X%*%(beta1[(2+n.fixed.pred):(n.pred+1+n.fixed.pred),]*beta1[(2+n.fixed.pred):(n.pred+1+n.fixed.pred),]))
      X<-cbind(1, fixed.pred, pred)
      V<-diag(rep(vy, times=N))+na.exclude(me.response)+ obs_var_con - diag(as.numeric(me.cov%*%(2*beta1[(2+n.fixed.pred):(n.pred+n.fixed.pred+1),]))) -diag(as.numeric(me.fixed.cov%*%(2*beta1[2:(length(beta1)-n.pred),])))
    }
    else
    {
      if(ultrametric==TRUE){
        s1<-as.numeric(s.X%*%(beta1[(2+n.fixed.pred):(n.pred+n.fixed.pred+1),]*beta1[(2+n.fixed.pred):(n.pred+n.fixed.pred+1),]))
        X<-cbind(1, fixed.pred, (1-(1-exp(-a*T))/(a*T))*pred)
        mcov<-diag(rowSums(matrix(data=as.numeric(me.cov)*t(kronecker(2*beta1[(2+n.fixed.pred):(n.pred+n.fixed.pred+1),], (1-(1-exp(-a*T))/(a*T)))), ncol=n.pred)))
      }
      else{
        s1<-as.numeric(s.X%*%(beta1[(4+n.fixed.pred):(n.pred+n.fixed.pred+3),]*beta1[(4+n.fixed.pred):(n.pred+n.fixed.pred+3),]))
        nu.X<-cbind(1-exp(-a*T), 1-exp(-a*T)-(1-(1-exp(-a*T))/(a*T)), exp(-a*T), fixed.pred, (1-(1-exp(-a*T))/(a*T))*pred)
        mcov<-diag(rowSums(matrix(data=as.numeric(me.cov)*t(kronecker(2*beta1[(4+n.fixed.pred):(n.pred+n.fixed.pred+3),], (1-(1-exp(-a*T))/(a*T)))), ncol=n.pred)))
        
      }

      

      cm1<-(s1/(2*a)+vy)*(1-exp(-2*a*ta))*exp(-a*tij)

      if(ultrametric==TRUE)
      {
        V<-cm1+(s1*ta*cm2)+na.exclude(me.response)+ obs_var_con-mcov -diag(as.numeric(me.fixed.cov%*%(2*beta1[2:(length(beta1)-n.pred),])))
      }
      else
      {
        V<-cm1+(s1*ta*cm2)+na.exclude(me.response)+ obs_var_con- mcov -diag(as.numeric(me.fixed.cov%*%(2*beta1[4:(length(beta1)-n.pred),])))
        
      }
    } # END OF ELSE CONDITION FOR HALF-LIFE = 0
    
    # INTERMEDIATE ESTIMATION OF OPTIMAL REGRESSION #
    
    V.inverse<-solve(V)
    if(hl==0)
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
  
  if(hl==0)
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
    sup1 <- -N/2*log(2*pi)-0.5*log(det(V))-0.5*(t(resid) %*% V.inverse%*%resid);
  }  # END OF CONDITION FOR HALF-LIFE = 0 #
  
  print(as.numeric(round(cbind(hl, vy, sup1, t(beta1)), 4)))
  
  list(support = sup1,
       V = V,
       beta1 = beta1,
       X = X,
       #beta1.var = beta.i.var,
       alpha.est = a,
       vy.est = vy)
}


   