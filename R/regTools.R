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
    if(is.null(fixed.fact)){
      if(!is.null(intercept)){
        "Intercept"
      }else{
        c("Intercept",
          if (!is.null(random.cov) | !is.null(fixed.cov))c("b0", "b1Xa") else NULL)
      }
    }else{
      if(is.null(intercept)){
        c("Ya", levels(as.factor(fixed.fact))[unique(factor(fixed.fact))])
      }else{
        levels(as.factor(fixed.fact))[unique(factor(fixed.fact))]
      }
    }
}