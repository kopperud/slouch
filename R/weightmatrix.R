parent <- function(phy, x){
  m <- which(phy$edge[, 2] == x)
  return(phy$edge[m, 1])
}

lineage.nodes <- function(phy, x){
  k <- x
  N <- length(phy$tip.label)
  while(x != N + 1){
    k <- c(k, parent(phy, x))
    x <- tail(k, n = 1)
  }
  return(k)
}

lineage.constructor <- function(phy, e, regimes){
  nodes <- lineage.nodes(phy, e)
  which.regimes = lapply(levels(regimes), function(x) {res <- match(regimes[nodes], x); res[is.na(res)] <- 0; return(res) })
  names(which.regimes) <- levels(regimes)
  
  nodes_time <-  ape::node.depth.edgelength(phy)[nodes]
  timeflip <- nodes_time[1] - nodes_time ## Time from tip to node(s)
  t_end <- tail(timeflip, n = -1) ## Time from tip to end of segment(s)
  t_beginning <- head(timeflip, n = -1) ## Time from tip to beginning of segment(s)
  
  return(list(nodes = nodes, 
              nodes_time = nodes_time,
              t_end = t_end,
              t_beginning = t_beginning,
              which.regimes = which.regimes))
}

weights_segments <- function(a, lineage){
  res <- c(exp(-a * lineage$t_beginning) - exp(-a * lineage$t_end), 
           exp(-a * lineage$nodes_time[1]))
  return(res)
}

weights_regimes <- function(a, lineage) {
  #nt <- lineage$nodes_time
  res <- weights_segments(a, lineage) ## Rcpp wrapper, equivalent to above commented out code
  w <- vapply(lineage$which.regimes, function(e) sum(e*res), FUN.VALUE = 0) ## Sum coefficients for which the respective regime is equal
  return(w)
}

weight.matrix <- function(phy, a, lineages){
  if(a > 300000000000) a <- 300000000000
  res <- t(vapply(lineages, function(x) weights_regimes(a, x), 
                  FUN.VALUE = numeric(length(lineages[[1]]$which.regimes))) ## Specify type of output
           )

  rownames(res) <- phy$tip.label
  return(res)
}

## Thanks to user "snaut" at stackoverflow, http://stackoverflow.com/users/1999873/snaut

concat.factor <- function(...){
  as.factor(do.call(c, lapply(list(...), as.character)))
}