#' Function to fit Ornstein-Uhlenbeck models of trait evolution
#'
#' 
#' @param phy an object of class 'phylo', must be rooted.
#' @param species a character vector of species tip labels, typically the "species" column in a data frame. This column needs to be an exact match and same order as phy$tip.label
#' @param sigma2_y_values alternative to vy_values, if the stationary variance is reparameterized as the variance parameter for the Brownian motion.
#' @param hl_values a vector of candidate phylogenetic half-life values to be evaluated in grid search. Optional.
#' @param vy_values a vector of candidate stationary variances for the response trait, to be evaluated in grid search. Optional.
#' @param estimate.Ya a logical value indicathing whether "Ya" should be estimated. If true, the intercept K = 1 is expanded to Ya = exp(-a*t) and b0 = 1-exp(-a*t). If models with categorical covariates are used, this will instead estimate a separate primary optimum for the root niche, "Ya". This only makes sense for non-ultrametric trees. If the tree is ultrametric, the model matrix becomes singular.
#' @param estimate.bXa a logical value indicathing whether "bXa" should be estimated. If true, bXa = 1-exp(-a*t) - (1-(1-exp(-a*t))/(a*t)) is added to the model matrix, estimating b*Xa. Same requirements as for estimating Ya.
#' @param response a numeric vector of a trait to be treated as response variable
#' @param mv.response numeric vector of the observational variances of each response trait. E.g if response is a mean trait value, mv.response is the within-species squared standard error of the mean.
#' @param fixed.fact factor of regimes on the terminal edges of the tree, in same order as species. If this is used, phy$node.label needs to be filled with the corresponding internal node regimes, in the order of node indices (root: n+1),(n+2),(n+3), ...
#' @param direct.cov Direct effect independent variables
#' @param mv.direct.cov Estimation variances for direct effect independent variables. Must be the same shape as direct.cov
#' @param mcov.direct.cov Estimation covariances between the response variable and direct effect independent variables. Most be the same shape as direct.cov
#' @param random.cov Independent variables each modeled as a brownian motion
#' @param mv.random.cov Estimation variances for the brownian covariates. Must be the same shape as random.cov
#' @param mcov.random.cov Estimation covariances between the response variable and random effect independent variables. Most be the same shape as random.cov
#' @param hessian use the approximate hessian matrix at the likelihood peak as found by the hillclimber, to compute standard errors for the parameters that enter in parameter search.
#' @param support a scalar indicating the size of the support set, defaults to 2 units of log-likelihood.
#' @param convergence threshold of iterative GLS estimation for when beta is considered to be converged.
#' @param nCores number of CPU cores used in grid-search. If 2 or more cores are used, all print statements are silenced during grid search. If performance is critical it is recommended to compile and link R to a multithreaded BLAS, since most of the heavy computations are common matrix operations. Even if a singlethreaded BLAS is used, this may or may not improve performance, and performance may vary with OS.
#' @param hillclimb logical, whether to use hillclimb parameter estimation routine or not. This routine (L-BFGS-B from optim()) may be combined with the grid-search, in which case it will on default start on the sigma and halflife for the local ML found by the grid-search.
#' @param lower lower bounds for the optimization routine, defaults to c(0,0). First entry in vector is half-life, second is stationary variance. When running direct effect models without observational error, it may be useful to specify a positive lower bounds for the stationary variance, e.g c(0, 0.001), since the residual variance-covariance matrix is degenerate when sigma = 0.
#' @param upper upper bounds for the optimization routine, defaults to c(Inf, Inf).
#' @param verbose a logical value indicating whether to print a summary in each iteration of parameter search. May be useful when diagnosing unexpected behaviour or crashes.
#' 
#' 
#' @return An object of class 'slouch', essentially a list with the following fields:
#' 
#' \item{parameter_space}{a list of the entire parameter space traversed by the grid search and the hillclimber as applicable.}
#' \item{tree}{a list of parameters concerning the tree:
#' \itemize{
#' \item{phy - an object of class 'phy'}
#' \item{T.term - a numeric vector including the time from the root of the tree to the tip, for all taxa 1,2,3... n.}
#' \item{ta - for all pairs of species, the time from their most recent common ancestor (mrca) to the root of the tree.}
#' \item{tia - for all pairs of species, the time from their mrca to the tip of species i.}
#' \item{tja - the transpose of tia.}
#' \item{tij - for all pairs of species, the time from species i to their mrca, plus the time from their mrca to species j. In other words, tia + transpose(tia).}
#' \item{times - for all nodes (1,2,3... n, root, root+1, ...) in the tree, the time from the root to said node.} 
#' \item{lineages - for all species (1,2,3... n), a list of their branch times and regimes as painted on the tree.}
#' \item{regimes - for all nodes (1,2,3... n, root, root+1, ...) in the tree, the respective regime as specified by "\code{phy$node.label}" and "\code{fixed.fact}".}
#' 
#' }
#' }
#' 
#' \item{modfit}{a list of statistics to characterize model fit}
#' \item{supportplot}{a list or matrix used to plot the grid search}
#' \item{supported_range}{a matrix indicating the interval of grid search that is within the support region. If the grid search values are carefully selected, this may be used to estimate the true support region.}
#' \item{V}{the residual variance-covariance matrix for the maximum likelihood model as found by parameter search.}
#' \item{evolpar}{maximum likelihood estimates of parameters under the chosen model.}
#' \item{beta_primary}{regression coefficients and associated objects. Whether the regression coefficients are to be interpreted as optima or not depend on the type of model and model estimates.}
#' \item{beta_evolutionary}{under a random effect model, "beta_evolutionary" is the evolutionary regression coefficients and associated objects.}
#' \item{n.par}{number of free parameters with which the likelihood criteria are penalized.}
#' \item{brownian_predictors}{under a random effect model, a matrix of means and standard errors for the independent Brownian motion variable(s). Not to be confused with the regression coefficients when the residuals are under a "bm" model.}
#' \item{climblog_df}{a matrix of the path trajectory of the hillclimber routine.}
#' \item{fixed.fact}{the respective regimes for all species (1,2,3... n).}
#' \item{control}{internal parameters for control flow.}
#' 
#' 
#' @export
slouch.fit <- function(phy,
                       species = NULL,
                       hl_values = NULL, 
                       vy_values = NULL, 
                       sigma2_y_values = NULL,
                       response, 
                       mv.response=NULL, 
                       fixed.fact=NULL,
                       direct.cov=NULL, 
                       mv.direct.cov=NULL, 
                       mcov.direct.cov=NULL, 
                       random.cov=NULL, 
                       mv.random.cov=NULL, 
                       mcov.random.cov=NULL,
                       estimate.Ya = FALSE,
                       estimate.bXa = FALSE,
                       hessian = F,
                       support = 2, 
                       convergence = 0.000001,
                       nCores = 1,
                       hillclimb = FALSE,
                       lower = c(1e-8, 1e-8),
                       upper = Inf,
                       verbose = FALSE)
{
  if((sum(c(is.null(hl_values), is.null(vy_values), is.null(sigma2_y_values))) > 1) & !hillclimb){
    stop("Choose at minimum a 1x1 grid, or use the hillclimber routine.")
  }
  if(!is.null(vy_values) & !is.null(sigma2_y_values)){
    stop("Choose either \"vy_values\" or \"sigma2_y_values\", not both.")
  }
  
  if(!is.null(random.cov)){
    if(ncol(as.matrix(random.cov))==1) {
      names.random.cov <- deparse(substitute(random.cov))
    }else{
      names.random.cov <- colnames(random.cov)
    }
  }
  
  if(!is.null(direct.cov)){
    if(ncol(as.matrix(direct.cov))==1) {
      names.direct.cov <- deparse(substitute(direct.cov))
    }else{
      names.direct.cov <- colnames(direct.cov)
    }
  }
  
  .slouch.fit(phy = phy,
              species = species,
              hl_values = hl_values, 
              vy_values = vy_values, 
              sigma2_y_values = sigma2_y_values,
              response = response, 
              mv.response = mv.response, 
              fixed.fact = fixed.fact,
              direct.cov = direct.cov, 
              mv.direct.cov = mv.direct.cov, 
              mcov.direct.cov = mcov.direct.cov, 
              random.cov = random.cov, 
              mv.random.cov = mv.random.cov, 
              mcov.random.cov = mcov.random.cov,
              estimate.Ya = estimate.Ya,
              estimate.bXa = estimate.bXa,
              hessian = hessian,
              model = "ou",
              support = support, 
              convergence = convergence,
              nCores = nCores,
              hillclimb = hillclimb,
              lower = lower,
              upper = upper,
              verbose = verbose,
              names.direct.cov = names.direct.cov,
              names.random.cov = names.random.cov)
}

#' Function to fit Brownian-motion models of trait evolution
#'
#' 
#' 
#' @inherit slouch.fit return
#'
#' @inheritParams slouch.fit
#' 
#' @param sigma2_y_values a vector of one or more candidates for sigma squared (y) to be evaluated in grid search.
#' @param estimate.Ya independently estimates the ancestral state under Brownian motion. Note that, for an intercept model, the intercept IS the ancestral state estimate (since there are no directional or stabilizing trends in a standard Brownian motion).
#' @param lower lower bounds for the optimization routine, defaults to 1e-8. When running direct effect models without observational error, it may be useful to specify a positive lower bounds for the sigma squared, since the residual variance-covariance matrix is degenerate when sigma = 0.
#' @param upper upper bounds for the optimization routine, defaults to 10 * var(response) * max(treeheight).
#' 
#' @export
brown.fit <- function(phy,
                      species = NULL,
                      sigma2_y_values = NULL,
                      response, 
                      mv.response = NULL, 
                      fixed.fact = NULL,
                      direct.cov = NULL, 
                      mv.direct.cov = NULL, 
                      mcov.direct.cov = NULL, 
                      random.cov = NULL, 
                      mv.random.cov = NULL, 
                      mcov.random.cov = NULL,
                      estimate.Ya = FALSE,
                      hessian = FALSE,
                      support = 2, 
                      convergence = 0.000001,
                      nCores = 1,
                      hillclimb = FALSE,
                      lower = 1e-8,
                      upper = NULL,
                      verbose = FALSE)
{
  if(is.null(sigma2_y_values) & !hillclimb){
    stop("Choose at minimum a sigma2_y_values of length one or use the hillclimber routine.")
  }
  if(is.null(upper)){
    upper <- 10 * stats::var(response) / max(ape::node.depth.edgelength(phy))
  }
  
  if(!is.null(random.cov)){
    if(ncol(as.matrix(random.cov))==1) {
      names.random.cov <- deparse(substitute(random.cov))
    }else{
      names.random.cov <- colnames(random.cov)
    }
  }
  
  if(!is.null(direct.cov)){
    if(ncol(as.matrix(direct.cov))==1) {
      names.direct.cov <- deparse(substitute(direct.cov))
    }else{
      names.direct.cov <- colnames(direct.cov)
    }
  }
  
  .slouch.fit(phy = phy,
              species = species,
              sigma2_y_values = sigma2_y_values,
              response = response, 
              mv.response = mv.response, 
              fixed.fact = fixed.fact,
              direct.cov = direct.cov, 
              mv.direct.cov = mv.direct.cov, 
              mcov.direct.cov = mcov.direct.cov, 
              random.cov = random.cov, 
              mv.random.cov = mv.random.cov, 
              mcov.random.cov = mcov.random.cov,
              estimate.Ya = estimate.Ya,
              estimate.bXa = FALSE,
              hessian = hessian,
              model = "bm",
              support = support, 
              convergence = convergence,
              nCores = nCores,
              hillclimb = hillclimb,
              lower = lower,
              upper = upper,
              verbose = verbose,
              names.direct.cov = names.direct.cov,
              names.random.cov = names.random.cov)
}


#' Logarithmically spaced sequence
#'
#' @param from the starting value of the sequence. Must be positive.
#' @param to the terminal value of the sequence. Must be larger than input to "from".
#' @param length.out desired length of the sequence. Must not be negative.
#'
#' @return A sequence of logarithmically spaced numbers.
#' @export
#'
#' @examples
#' lseq(1, 1000, length.out = 4)
lseq <- function(from=1, to=100000, length.out=6) {
  # logarithmic spaced sequence
  # blatantly stolen from library("emdbook"), because need only this
  exp(seq(log(from), log(to), length.out = length.out))
}