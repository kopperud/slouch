######## THIS FILE #########
## Extra helper functions ##

## Calculate determinant of V
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

# Find beta coef names
coef.names.factor <- function(modelpar, treepar, seed){
  names.intercept <- 
    if(is.null(modelpar$fixed.fact)){
      if(treepar$ultrametric | is.null(modelpar$random.cov)){
        "Intercept"
      }else{
        c("Ya", "Xa", "Bo")
      }
    }else{
      if(is.null(modelpar$intercept)){
        c("Ya", levels(modelpar$regime.specs)[unique(modelpar$regime.specs)])
      }
      levels(modelpar$regime.specs)[unique(modelpar$regime.specs)]
    } 
  # names.continuous <- c(if(seed$n.fixed.pred==1) deparse(substitute(modelpar$fixed.cov)) else colnames(modelpar$fixed.cov), 
  #                       if(seed$n.pred == 1) deparse(substitute(modelpar$random.cov)) else colnames(modelpar$random.cov))
  #c(names.intercept, names.continuous)
  names.intercept
}
