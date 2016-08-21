#### THIS FILE #######
#### Extra functions, mostly for regression
####
####

# Part "two" of variance-covariance matrix, hansen et al. 2008. Everything after sigma^2_theta * ta ????? Looks like equation is modified to 
# account for non-ultrametric trees.

make.cm2 <- function(a,tia,tja,ta,N,T.term){
  T.row <- replicate(N,T.term)
  T.col <- t(T.row)
  num.prob <- ifelse(ta == 0, 1, (1-exp(-a*ta))/(a*ta))
  #common_term <- 1-exp(-a*T.row)
  
  return(((1-exp(-a*T.row))/(a*T.row))*((1-exp(-a*T.col))/(a*T.col)) - (exp(-a*tia)*(1-exp(-a*T.row))/(a*T.col) + exp(-a*tja)*(1-exp(-a*T.row))/(a*T.row))*num.prob)
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

# ## General test for beta convergence
# test.conv <- function(beta.i, beta1, convergence, con.count, ultrametric){
#   if(ultrametric | (is.null(fixed.cov) & is.null(random.cov))){
#     test <- ifelse(abs(as.numeric(beta.i - beta1)) <= convergence, 0, 1)
#   }else{
#     test <- ifelse(abs(as.numeric(beta.i - beta1))[-(1:2)] <= convergence, 0, 1)
#     ## Effectively removes beta[1:2] <= 0.001 from being criteria in convergence, when non-ultrametric with continuous covariates
#   }
#   
#   if(sum(test)==0) return (TRUE)
#   if(con.count >= 50)
#   {
#     message("Warning, estimates did not converge after 50 iterations, last estimates printed out")
#     return(TRUE)
#   }
#   return(FALSE)
# }