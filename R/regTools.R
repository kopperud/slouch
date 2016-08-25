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

# Find beta coef names of intercepts + factor variables
coef.names.factor <- function(fixed.fact, random.cov, fixed.cov, intercept){
  names.intercept <- 
    if(is.null(fixed.fact)){
      #if(treepar$ultrametric | is.null(modelpar$random.cov)){
      if(!is.null(intercept)){
        "Intercept"
      }else{
        c("Intercept", "b0", if (!is.null(random.cov) | !is.null(fixed.cov))"b1Xa" else NULL)
      }
    }else{
      if(is.null(intercept)){
        c("Ya", levels(regime.specs)[unique(factor(fixed.fact))])
      }
    } 
  names.intercept
}