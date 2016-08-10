#### THIS FILE #######
#### Extra functions, mostly for regression
####
####

# Part "two" of variance-covariance matrix, hansen et al. 2008. Everything after sigma^2_theta * ta ????? Looks like equation is modified to 
# account for non-ultrametric trees. Confirm?

make.cm2 <- function(a,tia,tja,ta,N,T.term){
  T.row <- replicate(N,T.term)
  T.col <- t(T.row)
  num.prob <- ifelse(ta == 0, 1, (1-exp(-a*ta))/(a*ta))
  common_term <- 1-exp(-a*T.row)
  
  return(((common_term)/(a*T.row))*((1-exp(-a*T.col))/(a*T.col)) - (exp(-a*tia)*(common_term)/(a*T.col) + exp(-a*tja)*(common_term)/(a*T.row))*num.prob)
}

test.conv.rReg <- function(beta.i, beta1, n.pred, convergence, con.count, ultrametric){
  if (ultrametric == TRUE) {
    fstart <- 0
    y <- 1
  }
  else {
    fstart <- 3
    y <- 3
  }
  test<-matrix(nrow=(n.pred + y))
  for(f in (1 + fstart):(n.pred + y))
  {
    if(abs(as.numeric(beta.i[f]-beta1[f]))<=convergence){
      test[(f - fstart)]=0
    }else {
      test[(f - fstart)]=1
    }
  }
  print(beta1); print(beta.i); print(n.pred)
  if(sum(test)==0) return (TRUE)
  if(con.count >= 50)
  {
    message("Warning, estimates did not converge after 50 iterations, last estimates printed out")
    return(TRUE)
  }
  return(FALSE)
}

test.conv.fReg <- function(beta.i = beta.i, beta1 = beta1, convergence = convergence, con.count = con.count, ultrametric = ultrametric){
  test<-matrix(nrow=length(beta.i))
  for(f in 1:(length(beta.i)))
  {
    if(abs(as.numeric(beta.i[f]-beta1[f]))<=convergence){
      test[f]=0
    }else {
      test[f]=1
    }
  }
  if(sum(test)==0) return (TRUE)
  if(con.count >= 50)
  {
    message("Warning, estimates did not converge after 50 iterations, last estimates printed out")
    return(TRUE)
  }
  return(FALSE)
}

mk.log.det.V <- function(V, N){
  det.V<-det(V)
  if(det.V==0){
    print(paste("Warning: Determinant of V = 0"))
    #Minimum value of diagonal scaling factor
    inv.min.diag.V<-1/min(diag(V))
    V<-V*inv.min.diag.V
    #Rescale and log determinant
    log(det(V))+log(min(diag(V)))*N
  }
  else {
    log(det.V)
  }
}

mk.obs_var_con <- function(a, hl, beta1, T, N, xx, x.ols, error_condition){
  if (hl == 0){
    y <- 1
  }   else {
    y <- ((1-(1-exp(-a*T))/(a*T))*(1-(1-exp(-a*T))/(a*T)))
  }
  obs_var_con <- matrix(0, nrow=N, ncol=N)
  for (e in seq(from=1, to=ncol(x.ols), by=1)){
    for (j in seq(from=1, to=ncol(x.ols), by=1)) {
      tmp <- error_condition[xx[e]:(e*N),xx[j]:(j*N)]*beta1[e]*beta1[j]*y
      obs_var_con <- obs_var_con + tmp
    }
  }
  obs_var_con
}



# Deprecated.
# make.beta1.rReg <- function(hl, x.ols, Y, ultrametric){
#   if (hl != 0 & ultrametric == FALSE){
#     rbind(0, 0, solve(t(x.ols)%*%x.ols)%*%(t(x.ols)%*%Y))
#   } else{
#     solve(t(x.ols)%*%x.ols)%*%(t(x.ols)%*%Y)
#   }
# }

