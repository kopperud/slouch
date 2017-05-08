## Methods




#' Print slouch-objects
#' 
#' @param x An object of class 'slouch'
#' @param ... Additional arguments, unused.
#' 
#' @export
print.slouch <- function(x, ...){
  message("Important - Always inspect the likelihood surface of the model parameters in the 3D-grid
          before evaluating model fit & results. If the likelihood search space does not contain
          the true maximum likelihood, the model outputs will reflect this.")
  message("")
  message("Model parameters")
  print(x$oupar)
  
  if (!is.null(x$hlvy_grid_interval)){
    message("Interval of parameters in 3d plot (Very sensitive to grid mesh, grid size and local ML estimate)")
    print(x$hlvy_grid_interval)
  }

  message("Optimal regression")
  print(x$opt.reg$coefficients)

  message("Optimal regression - bias-corrected")
  print(x$opt.reg$coefficients_bias_corr)

  if (!is.null(x$ev.reg)){
    message("Evolutionary regression")
    print(x$ev.reg$coefficients)
  }
  
  if (!is.null(x$brownian_predictors)){
    message("Stochastic predictor")
    print(x$brownian_predictors)
  }
  
  message("Model fit")
  print(x$modfit)
  
  plot(x)
}


#' Plot slouch-objects
#' @description Graphical plot of parameter space traversed in order to find ML-estimate of the model.
#'
#' @param x An object of class 'slouch'
#' @param ... Additional parameters passed to persp()
#'
#' @export
plot.slouch <- function(x, ...){
  plotpars <- graphics::par()
  if (!is.null(x$supportplot) & !is.null(x$climblog_matrix)){
    graphics::par(mfrow = c(1, 2))
  }

  if (!is.null(x$supportplot)){
    graphics::persp(x$supportplot$hl,
          x$supportplot$vy,
          x$supportplot$z,
          theta = 30,
          phi = 30, 
          expand = 0.5, 
          col = "NA",
          ltheta = 120, 
          shade = 0.75,
          main = "Grid search",
          ticktype = "detailed",
          xlab = "Phylogenetic half-life", 
          ylab = "Stationary variance", 
          zlab = "Log-likelihood", 
          zlim = c(min(x$supportplot$z), 0),
          ...)
  }
  
  if (!is.null(x$climblog_matrix)){
    hl <- x[["climblog_matrix"]][["hl"]]
    vy <- x[["climblog_matrix"]][["vy"]]
    index <- x[["climblog_matrix"]][["index"]]
    
    graphics::plot(x = hl, 
         y = vy, 
         main = "Path of hillclimber",
         col = grDevices::gray.colors(length(index),
                                      start = 0.8, end = 0.05, gamma = 1)[index],
         pch = 19,
         ylim = c(min(vy), max(vy)),
         xlim = c(min(hl), max(hl)),
         xlab = "Phylogenetic half-life",
         ylab = "Stationary variance")
    text(hl[1], vy[1], "Start")
    text(hl[length(hl)], vy[length(vy)], "End")
  }

  graphics::par(mfrow = plotpars$mfrow)
}


#' Extract Log-Likelihood
#'
#' @param object An object of class 'slouch'
#' @param ... Additional arguments.
#'
#' @return An object of class 'logLik'
#' @export
logLik.slouch <- function(object, ...){
  return(structure(object$modfit[1, 1], df = object$n.par, class = "logLik"))
}




#' Plot the internal regimes for a given fitted model
#'
#' @param x an object of class 'slouch'
#' @param ... additional parameters
#'
#' @return nothing
#' @export
regimeplot.slouch <- function(x, ...){
  stopifnot(!is.null(x$fixed.fact))
  
  regimes <- c(x$fixed.fact, x$tree$phy$node.label)
  regimes_edges <- regimes[x$tree$phy$edge[,2]]
  
  plot(x$tree$phy, edge.col = factor(regimes_edges))
}

#' Title
#'
#' @param x an object of class 'slouch'
#' @param ... additional parameters
#'
#' @return
#' @export
regimeplot <- function(x, ...){
  UseMethod("regimeplot")
}