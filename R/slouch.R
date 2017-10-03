#'SLOUCH: Stochastic Linear Ornstein Uhlenbeck Comparative Hypotheses
#'
#'
#'
#'
#'@importFrom stats lm.fit na.exclude optim runif var logLik
#'@importFrom utils tail head
#'@importFrom grDevices palette gray.colors
#'@importFrom graphics plot text points axis
#'
#' @section References:
#'
#' \itemize{
#'  \item Hansen, T. F. (1997). Stabilizing Selection and the Comparative Analysis of Adaptation. Evolution, 51(5), 1341. https://doi.org/10.2307/2411186
#'  
#'  \item Hansen, T. F., Pienaar, J., & Orzack, S. H. (2008). A comparative method for studying adaptation to a randomly evolving environment. Evolution, 62(8), 1965–1977. https://doi.org/10.1111/j.1558-5646.2008.00412.x
#'  
#'  \item Labra, A., Pienaar, J., & Hansen, T. F. (2009). Evolution of Thermal Physiology in Liolaemus Lizards: Adaptation, Phylogenetic Inertia, and Niche Tracking. The American Naturalist, 174(2), 204–220. https://doi.org/10.1086/600088
#'  
#'  \item Hansen, T. F., & Bartoszek, K. (2012). Interpreting the evolutionary regression: The interplay between observational and biological errors in phylogenetic comparative studies. Systematic Biology, 61(3), 413–425. https://doi.org/10.1093/sysbio/syr122
#'  
#'  \item Escudero, M., Hipp, A. L., Hansen, T. F., Voje, K. L., & Luceño, M. (2012). Selection and inertia in the evolution of holocentric chromosomes in sedges (Carex, Cyperaceae). New Phytologist, 195(1), 237–247. https://doi.org/10.1111/j.1469-8137.2012.04137.x
#' }
"_PACKAGE"