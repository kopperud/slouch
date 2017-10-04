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
  if(length(x$evolpar) > 1){
    evolpar <- matrix(x$evolpar, ncol = 1, dimnames=list(c("Rate of adaptation", "Phylogenetic half-life", "Stationary variance", "Sigma squared (y)","Phylogenetic correction factor"), "Estimate"))
  }else{
    evolpar <- matrix(x$evolpar, dimnames = list(c("Brownian sigma2 (y)"), "Estimate"))
  }
  print(evolpar)
  
  if (!is.null(x$supported_range)){
    message("Interval of parameters in 3d plot (Very sensitive to grid mesh, grid size and local ML estimate)")
    print(x$supported_range)
  }

  if(x$control$model == "ou") message("Optimal regression") else message("Phylogenetic regression")
  print(x$opt.reg$coefficients)

  if(x$control$model == "ou") message("Optimal regression - bias-corrected") else message("Phylogenetic regression - bias-corrected")
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
  m <- as.matrix(x$modfit)
  colnames(m) <- "Values"
  print(m)
  
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
  if (!is.null(x$supportplot) & !is.null(x$climblog_df)){
    graphics::par(mfrow = c(1, 2))
  }

  
  if (!is.null(x$supportplot)){
    if(x$control$model == "ou"){
      s3d <- if(!is.null(x$supportplot$vy)) x$supportplot$vy else x$supportplot$sigma2_y
      graphics::persp(x$supportplot$hl,
                      s3d,
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
                      ylab = if(!is.null(x$supportplot$vy)) "Stationary variance" else "Sigma squared (y)", 
                      zlab = "Log-likelihood", 
                      zlim = c(min(x$supportplot$z), 0),
                      ...)
    }else{
      hline <- logLik(x) - x$control$support
      graphics::plot(x = log(x$supportplot$sigma2_y),
                     y = x$supportplot$loglik,
                     pch = 19,
                     xlab = "sigma squared (y)",
                     ylab = "logL",
                     xaxt = "n",
                     ylim = c(min(hline, min(x$supportplot$loglik)), 
                              max(c(x$supportplot$loglik, logLik(x)))),
                     main = "Grid search")
      graphics::axis(1, 
                     at = log(x$supportplot$sigma2_y), 
                     labels = format(x$supportplot$sigma2_y, digits = 1),
                     las = 1)
      
      graphics::abline(h = hline, lwd = 2, col = "red")
    }
  }
  
  if (!is.null(x$climblog_df)){
    logL <- x$climblog_df$loglik
    hl <- x$climblog_df$hl
    vy <- x$climblog_df$vy
    sigma2_y <- x$climblog_df$sigma2_y
    index <- x$climblog_df$index
    s <- if(!is.null(vy)) vy else sigma2_y
    
    if(x$control$model == "ou"){
      graphics::plot(x = hl, 
                     y = s, 
                     main = "Path of hillclimber",
                     col = grDevices::gray.colors(length(index),
                                                  start = 0.8, end = 0.05, gamma = 1)[index],
                     pch = 19,
                     ylim = c(min(s), max(s)),
                     xlim = c(min(hl), max(hl)),
                     xlab = "Phylogenetic half-life",
                     ylab = if(!is.null(vy)) "Stationary variance" else "Sigma squared (y)"
                     )
      
      points(hl[length(hl)], s[length(s)], col = "red", pch = 19)
      text(hl[1], s[1], "Start")
      text(hl[length(hl)], s[length(s)], "End")
    }else{
      graphics::plot(x = log(sigma2_y),
                     y = logL,
                     xlab = "sigma squared (y)", pch = 19,
                     xaxt = "n",
                     main = "Path of hillclimber")
      graphics::axis(1, 
                     at = log(x$climblog_df$sigma2_y), 
                     labels = round(x$climblog_df$sigma2_y, 2))
      points(log(tail(sigma2_y, n = 1)), 
             tail(logL, n = 1),
             pch = 19, col = "red")
      points(log(head(sigma2_y, n = 1)), 
             head(logL, n = 1),
             pch = 19, col = "darkblue")
      text(log(sigma2_y[1]), logL[1], "Start")
      text(tail(log(sigma2_y), n = 1), tail(logL, n = 1), "End")
    }
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
  return(structure(object$modfit$Support, df = object$n.par, class = "logLik"))
}




#' Plot the internal regimes for a given fitted model
#'
#' @param x an object of class 'slouch'
#' @param ... additional parameters passed to plot.phylo(...)
#'
#' @return nothing
#' @export
regimeplot.slouch <- function(x, ...){
  stopifnot(!is.null(x$fixed.fact))
  
  regimes <- concat.factor(x$fixed.fact, x$tree$phy$node.label)
  regimes_edges <- regimes[x$tree$phy$edge[,2]]
  
  ape::plot.phylo(x$tree$phy, edge.col = factor(regimes_edges), ...)
  
  p <- grDevices::palette()
  print("Colors:")
  for (e in seq_along(levels(regimes))){
    print(paste0(levels(regimes)[e], ": ", p[e]))
  }
}

#' Title
#'
#' @param x an object of class 'slouch'
#' @param ... additional parameters
#'
#' @return nothing
#' @export
regimeplot <- function(x, ...){
  UseMethod("regimeplot")
}