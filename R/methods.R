## Methods




#' Print, minimalist output
#' 
#' @param x An object of class 'slouch'
#' @param ... additional parameters
#'
#' @examples 
#' data(artiodactyla)
#' data(neocortex)
#' 
#' neocortex <- neocortex[match(artiodactyla$tip.label, neocortex$species), ]
#' 
#' m0 <- slouch.fit(phy = artiodactyla,
#'                  hl_values = seq(0.001, 4, length.out = 10),
#'                  vy_values = seq(0.001, 0.05, length.out = 10),
#'                  species = neocortex$species,
#'                  response = neocortex$neocortex_area_mm2_log_mean,
#'                  mv.response = neocortex$neocortex_se_squared,
#'                  random.cov = neocortex$brain_mass_g_log_mean,
#'                  mv.random.cov = neocortex$brain_se_squared,
#'                  fixed.fact = neocortex$diet,
#'                  hillclimb = FALSE)
#' @export
print.slouch <- function(x, ...){
  digits <- max(3, getOption("digits")- 3)
    
  cat(paste0(bold("Response: "), x$control$name.response, "\n"))
  cat("\n  ")
  
  cat(bold("Model fit:\n"))
  fit <- unlist(x$modfit[c("AICc", "Support", "R squared")])
  print(fit, digits = digits)
  
  cat("\n  ")
  cat(bold("ML estimate(s): \n"))
  
  evolpar <- unlist(x$evolpar)
  names(evolpar) <- parnames(names(evolpar))
  for (i in seq_along(evolpar)){
    cat(paste0(names(evolpar)[i], ":\t ", round(evolpar[i], digits), "\n"))
  }
  
  regpar <- x$beta_primary$coefficients[,1]
  names(regpar) <- rownames(x$beta_primary$coefficients)
  cat("\n  ")
  cat(bold("Coefficients: \n"))
  print(regpar, digits = digits)
}

#' Model Summary
#' 
#' @param object An object of class 'slouch'
#' @param ... Additional arguments, unused.
#'
#' @examples 
#' data(artiodactyla)
#' data(neocortex)
#' 
#' neocortex <- neocortex[match(artiodactyla$tip.label, neocortex$species), ]
#' 
#' m0 <- slouch.fit(phy = artiodactyla,
#'                  hl_values = seq(0.001, 4, length.out = 10),
#'                  vy_values = seq(0.001, 0.05, length.out = 10),
#'                  species = neocortex$species,
#'                  response = neocortex$neocortex_area_mm2_log_mean,
#'                  mv.response = neocortex$neocortex_se_squared,
#'                  random.cov = neocortex$brain_mass_g_log_mean,
#'                  mv.random.cov = neocortex$brain_se_squared,
#'                  fixed.fact = neocortex$diet,
#'                  hillclimb = FALSE)
#'                  
#' summary(m0) 
#' 
#' plot(m0, theta = 150)
#' @export
summary.slouch <- function(object, ...){
  x <- object
  cat(bold("Important"), "- Always inspect the likelihood surface of the model parameters with
          grid search before evaluating model fit & results.\n")
  cat("\n  ")
  cat(bold("Maximum-likelihood estimates\n"))
  evolpar <- as.matrix(x$evolpar, ncol = 1)
  rownames(evolpar) <- parnames(rownames(evolpar))
  if(!is.null(x$hessian)){
    evolpar <- cbind(evolpar, sqrt(diag(solve(-x$hessian))))
  }
  colnames(evolpar) <- c("Estimate", "Approximate Std. error (using hessian)")[1:ncol(evolpar)]
  
  if(!is.null(x$supportplot)){
    evolpar_names = c()
    for (m in names(x$evolpar)){
      par_range <- range(sapply(x$parameter_space, function(e) if(e$gridsearch) e$par[[m]] else NA), na.rm = T)
      
      p <- x$evolpar[[m]]
      
      if(p < par_range[1] | p > par_range[2]){
        evolpar_names = append(evolpar_names, paste(parnames(m), "(warning: outside grid)"))
      }else{
        evolpar_names = append(evolpar_names, parnames(m))
      }
    }
    rownames(evolpar) <- evolpar_names
  }
  
  print(evolpar)
  
  if (!is.null(x$brownian_predictors)){
    cat("\n  ")
    cat(bold("Stochastic predictor(s)\n"))
    print(x$brownian_predictors)
  }

  if(x$control$model == "ou"){
    if(is.null(x$evolpar$hl)){
      a <- x$evolpar$a
      hl <- log(2) / a
    }else{
      hl <- x$evolpar$hl
      a <- log(2) / hl      
    }
    
    if(!is.null(x$evolpar$vy)){
      vy <- x$evolpar$vy
      sigma2_y <- vy*2*a
    }else{
      vy <- x$evolpar$sigma2_y / (2 * a)
      sigma2_y <- x$evolpar$sigma2_y
    }
    
    rho <- mean(1 - (1 - exp(-a * x$tree$T.term)) / (a * x$tree$T.term))
    inferred <- list("Mean phylogenetic correction factor" = rho)
    
    if(!is.null(x$evolpar$hl)){
     inferred["Rate of adaptation"]  = a
    }
    if(!is.null(x$evolpar$a)){
      inferred["Phylogenetic half-life"]  = hl
    }
    
    if(!is.null(x$evolpar$vy)){
      inferred["Diffusion variance"] <- sigma2_y
    }
    if(!is.null(x$evolpar$sigma2_y)){
      inferred["Stationary variance"] <- vy
    }
    
  }else{
    inferred <- list()
    
    if(length(c(x$w_beta$which.regimes, x$w_beta$which.random.cov)) > 0){
      rho <- mean(x$tree$T.term / 2)
      inferred[["Mean phylogenetic correction factor"]] <- rho
    }
  }
  
  if (length(inferred) > 0){
    cat("\n  ")
    cat(bold("Inferred maximum-likelihood parameters\n"))
    inferred_matrix <- matrix(inferred, ncol = 1, dimnames = list(names(inferred), "Value"))
    print(inferred_matrix)
  }
  
  if (!is.null(x$supported_range)){
    cat("\n  ")
    cat(bold("Interval of parameters in 3d plot (Sensitive to grid mesh, grid size and local ML estimate)\n"))
    rownames(x$supported_range) <- parnames(rownames(x$supported_range))
    print(x$supported_range)
  }

  
  foo <- function(w, title){
    if(length(w) > 0){
      cat("\n  ")
      if(title == "Regime-dependent trends"){
        if(!x$control$estimate.Ya){
          title <- paste(title, "(assuming Ya = 0)")
        }
      }
      cat(bold(title));cat("\n")
      print(x$beta_primary$coefficients[w, ,drop=FALSE])
      
      if(title == "Optimal regression slope" | title == "Trend covariates"){
        cat("\n  ")
        cat(bold("Evolutionary regression slope\n"))
        print(x$beta_evolutionary$coefficients[w, ,drop = FALSE])
      }
      
      if(grepl("Regime trends", title)){
        cat("\n  ")
        cat(bold("Pairwise contrasts among trends:\n"))
        print(x$beta_primary$trend_diff)
      }
    }
  }
  
  if(x$control$model == "ou"){
    printnames <- c("Intercepts", "Regime optima", "Optimal regression slope", "Direct-effect slope")
  }else{
    printnames <- c("Intercepts", "Regime trends", "Trend covariates", "Direct-effect slope")
  }
  
  mapply(foo, w = x$w_beta, title = printnames)
  
  if(!is.null(x$beta_primary$K)){
    if(!is.diag(x$beta_primary$K)){
      cat("\n  ")
      cat(bold("Attenuation factor. Linear model coefficients (above) are not corrected for bias.\n"))
      m <- signif(x$beta_primary$K, 3)
      prmatrix(m, collab = rep_len("", ncol(m)))
    }
  }

  cat("\n  ")
  cat(bold("Model fit summary\n"))
  m <- as.matrix(sapply(x$modfit, signif, 3))
  colnames(m) <- "Values"
  print(m)
  
  plot(x)
}


#' Plot the hillclimber trajectory
#'
#' @param x An object of class 'slouch'
#' @param ... Additional arguments passed to 'plot.default(...)'
#'
#' @export
#' @examples 
#'library(slouch)
#'library(ape)
#'
#'data(neocortex)
#'data(artiodactyla)
#'
#'neocortex <- neocortex[match(artiodactyla$tip.label, neocortex$species), ]
#'
#'m0 <- slouch.fit(phy = artiodactyla,
#'                 species = neocortex$species,
#'                 response = neocortex$neocortex_area_mm2_log_mean,
#'                 mv.response = neocortex$neocortex_se_squared,
#'                 hillclimb = TRUE)
#'                 
#'hillclimbplot(m0)
#'
#'m1 <- brown.fit(phy = artiodactyla,
#'                species = neocortex$species,
#'                response = neocortex$neocortex_area_mm2_log_mean,
#'                mv.response = neocortex$neocortex_se_squared,
#'                hillclimb = TRUE)
#'                
#'hillclimbplot(m1)
hillclimbplot <- function (x, ...) {
  UseMethod("hillclimbplot", x)
}


#' @describeIn hillclimbplot Hillclimbplot for the 'slouch object'
#' @export
hillclimbplot.slouch <- function(x,...){
  
  if (!is.null(x$climblog_df)){
    logL <- x$climblog_df$loglik
    hl <- x$climblog_df$hl
    a <- x$climblog_df$a
    vy <- x$climblog_df$vy
    sigma2_y <- x$climblog_df$sigma2_y
    index <- x$climblog_df$index
    y_axis <- if(!is.null(vy)) vy else sigma2_y
    x_axis <- if(!is.null(a)) a else hl
    
    if(x$control$model == "ou"){
      graphics::plot(x = x_axis, 
                     y = y_axis, 
                     main = "Path of hillclimber",
                     col = grDevices::gray.colors(length(index),
                                                  start = 0.8, end = 0.05, gamma = 1)[index],
                     pch = 19,
                     ylim = range(y_axis),
                     xlim = range(x_axis),
                     xlab = if(!is.null(hl)) "Phylogenetic half-life" else "Rate of adaptation",
                     ylab = if(!is.null(vy)) "Stationary variance" else "Diffusion variance"
      )
      
      points(x_axis[length(x_axis)], y_axis[length(y_axis)], col = "red", pch = 19)
      text(x_axis[1], y_axis[1], "Start")
      text(x_axis[length(x_axis)], y_axis[length(y_axis)], "End")
      
    }else{
      graphics::plot(x = log(sigma2_y),
                     y = logL,
                     xlab = "Diffusion variance", pch = 19,
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
  }else{
    stop("Your model was not fitted using the hillclimber.")
  }
  
}


#' Plot Grid Search
#' @description Graphical plot of parameter space traversed by the grid search.
#'
#' @param x An object of class 'slouch'
#' @param ... Additional parameters passed to persp(...) or plot(...)
#' 
#' @inheritParams graphics::persp
#'
#' @examples 
#' 
#' data(artiodactyla)
#' data(neocortex)
#' 
#' neocortex <- neocortex[match(artiodactyla$tip.label, neocortex$species), ]
#' 
#' m0 <- slouch.fit(phy = artiodactyla,
#'                  hl_values = seq(0.001, 50, length.out = 15),
#'                  vy_values = seq(0.001, 3, length.out = 15),
#'                  species = neocortex$species,
#'                  response = neocortex$body_mass_g_log_mean,
#'                  mv.response = neocortex$body_mass_g_log_varmean,
#'                  fixed.fact = neocortex$diet)
#'                  
#' plot(m0)
#' @export
plot.slouch <- function(x, 
                        theta = 30,
                        phi = 30, 
                        expand = 0.5, 
                        shade = 0.75,
                        ...){

  if (!is.null(x$supportplot)){
    if(x$control$model == "ou"){
      s3d <- if(!is.null(x$supportplot$vy)) x$supportplot$vy else x$supportplot$sigma2_y
      x_axis <- if(!is.null(x$supportplot$hl)) x$supportplot$hl else x$supportplot$a

      if(all(dim(x$supportplot$z) > 1)){
        graphics::persp(x = x_axis,
                        y = s3d,
                        z = x$supportplot$z,
                        theta = theta,
                        phi = phi, 
                        expand = expand, 
                        ltheta = 120, 
                        shade = shade,
                        col = "NA",
                        main = "Grid search",
                        ticktype = "detailed",
                        xlab = if(!is.null(x$supportplot$hl)) "Phylogenetic half-life" else "Rate of adaptation", 
                        ylab = if(!is.null(x$supportplot$vy)) "Stationary variance" else "Diffusion variance", 
                        zlab = "Log-likelihood", 
                        zlim = c(min(x$supportplot$z), 0),
                        ...)
      }else{
        hline <- logLik(x) - x$control$support
        parnames <- c("hl", "a", "vy", "sigma2_y")
        which_one <- sapply(x$supportplot[parnames], function(e) length(e) == 1)
        which_variable <- sapply(x$supportplot[parnames], function(e) length(e) > 1)
        map <- stats::setNames(c("Phylogenetic half-life", "Rate of adaptation", "Stationary variance", "Diffusion variance"), parnames)
        
        xvar <- x$supportplot[parnames][which_variable]
        
        hline <- logLik(x) - x$control$support
        graphics::plot(x = xvar[[1]],
                       y = stats::na.omit(sapply(x$parameter_space, function(e) if(e$gridsearch) e$support else NA)),
                       pch = 19,
                       xlab = map[names(xvar)],
                       ylab = "logL",
                       xaxt = "n",
                       ylim = c(min(hline, min(x$supportplot$z)), 
                                max(c(x$supportplot$z, logLik(x)))),
                       main = paste0("Grid search", 
                                     " (at ", 
                                     parnames[which_one], 
                                     " = ", 
                                     x$supportplot[parnames][which_one][[1]], 
                                     ")"))
        graphics::axis(1, 
                       at = xvar[[1]], 
                       labels = format(xvar[[1]], digits = 1),
                       las = 1)
        graphics::abline(h = hline, lwd = 2, col = "red")
      }

    }else{
      hline <- logLik(x) - x$control$support
      graphics::plot(x = log(x$supportplot$sigma2_y),
                     y = x$supportplot$loglik,
                     pch = 19,
                     xlab = "Diffusion variance",
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
  
}


#' Extract Log-Likelihood
#'
#' @param object An object of class 'slouch'
#' @param ... Additional arguments.
#'
#' @return An object of class 'logLik'
#' @examples 
#' data(artiodactyla)
#' data(neocortex)
#' 
#' neocortex <- neocortex[match(artiodactyla$tip.label, neocortex$species), ]
#' 
#' m0 <- slouch.fit(phy = artiodactyla,
#'                  species = neocortex$species,
#'                  response = neocortex$body_mass_g_log_mean,
#'                  mv.response = neocortex$body_mass_g_log_varmean,
#'                  fixed.fact = neocortex$diet,
#'                  hillclimb = TRUE)
#'                  
#' logLik(m0)
#' @export
logLik.slouch <- function(object, ...){
  return(structure(object$modfit$Support, df = object$n.par, class = "logLik"))
}




#' @describeIn regimeplot Regimeplot for the 'slouch' object
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

#' Plot the internal regimes for a given fitted model
#'
#' @param x an object of class 'slouch'
#' @param ... additional parameters passed to plot.phylo(...)
#'
#' @return nothing
#' @examples 
#' 
#' data(artiodactyla)
#' data(neocortex)
#' 
#' neocortex <- neocortex[match(artiodactyla$tip.label, neocortex$species), ]
#' 
#' m0 <- slouch.fit(phy = artiodactyla,
#'                  species = neocortex$species,
#'                  response = neocortex$body_mass_g_log_mean,
#'                  mv.response = neocortex$body_mass_g_log_varmean,
#'                  fixed.fact = neocortex$diet,
#'                  hillclimb = TRUE)
#'                  
#' regimeplot(m0)
#' 
#' @export
regimeplot <- function(x, ...){
  UseMethod("regimeplot")
}


parnames <- function(x){
  x[x == "vy"] <- "Stationary variance"
  x[x == "sigma2_y"] <- "Diffusion variance"
  x[x == "hl"] <- "Phylogenetic half-life"
  x[x == "a"] <- "Rate of adaptation"
  return(x)
}