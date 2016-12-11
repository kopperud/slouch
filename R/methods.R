## Methods

#' Title
#'
#' @param x 
#'
#' @return
#' @export
#'
#' @examples
print.slouch <- function(x){
  message("Model parameters")
  print(x$oupar)
  
  if(!is.null(x$hlvy_grid_interval)){
    message("Interval of parameters in 3d plot")
    print(x$hlvy_grid_interval)
  }
  
  message("Optimal regression")
  print(x$opt.reg)
  
  if(!is.null(x$ev.reg)){
    message("Evolutionary regression")
    print(x$ev.reg)
  }
  
  if(!is.null(x$brownian_predictors)){
    message("Stochastic predictor")
    print(x$brownian_predictors)
  }
  
  message("Model fit")
  print(x$modfit)
  
  if(!is.null(x$supportplot)){
    plot(x)
  }
}

#' Title
#'
#' @param e 
#'
#' @return
#' @export
#'
#' @examples
plot.slouch <- function(e, theta = 30, ...){
  if (is.null(e$supportplot)){
    stop("Support grid not included.")
  }
  
  x <- e[["supportplot"]][["x"]]
  y <- e[["supportplot"]][["y"]]
  z <- e[["supportplot"]][["z"]]
  
  persp(x, y, z, theta = theta, phi = 30, expand = 0.5, col = "NA",
        ltheta = 120, shade = 0.75, ticktype = "detailed",
        xlab = "Phylogenetic half-life", ylab = "Stationary variance", zlab = "Log-likelihood", ...)
}