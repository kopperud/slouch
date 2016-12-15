## Methods




#' Title
#'
#' @param x 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
print.slouch <- function(x, ...){
  message("Important - Always inspect the likelihood surface of the model parameters in the 3D-grid before evaluating model fit & results. If the likelihood search space does not contain the true maximum likelihood, the model outputs will reflect this.")
  message("")
  message("Model parameters")
  print(x$oupar)
  
  if(!is.null(x$hlvy_grid_interval)){
    message("Interval of parameters in 3d plot")
    print(x$hlvy_grid_interval)
  }
  
  message("Optimal regression")
  print(x$opt.reg$coefficients)
  
  message("Optimal regression - bias-corrected")
  print(x$opt.reg$coefficients_bias_corr)
  
  if(!is.null(x$ev.reg)){
    message("Evolutionary regression")
    print(x$ev.reg$coefficients)
  }
  
  if(!is.null(x$brownian_predictors)){
    message("Stochastic predictor")
    print(x$brownian_predictors)
  }
  
  message("Model fit")
  print(x$modfit)
  
  plot(x)
}


#' Title
#'
#' @param x 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
plot.slouch <- function(x, ...){
  if (!is.null(x$supportplot)){
    #stop("Support grid not included.")
    x1 <- x[["supportplot"]][["x"]]
    y1 <- x[["supportplot"]][["y"]]
    z1 <- x[["supportplot"]][["z"]]
    
    persp(x1, y1, z1, theta = 30, phi = 30, expand = 0.5, col = "NA",
          ltheta = 120, shade = 0.75, ticktype = "detailed",
          xlab = "Phylogenetic half-life", ylab = "Stationary variance", zlab = "Log-likelihood", ...)
  }
  

  
  if (!is.null(x$climblog_matrix)){
    hl <- x[["climblog_matrix"]][["hl"]]
    vy <- x[["climblog_matrix"]][["vy"]]
    index <- x[["climblog_matrix"]][["index"]]
    
    plot(x = hl, 
         y = vy, 
         main ="Path of hillclimber",
         col = grDevices::gray.colors(length(index),
                                      start = 0.8, end=0.05, gamma = 1)[index],
         pch = 19,
         ylim = c(0, max(vy)),
         xlim = c(0, max(hl)),
         xlab = "Phylogenetic half-life",
         ylab = "Stationary variance")
    text(hl[1], vy[1], "Start")
    text(hl[length(hl)], vy[length(vy)], "End", ...)
  }

}


#' Title
#'
#' @param object 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
logLik.slouch <- function(object, ...){
  return(structure(object$modfit[1,1], df = object$n.par, class = "logLik"))
}