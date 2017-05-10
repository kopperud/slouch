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
  
  return(list(nodes = nodes, 
              nodes_time = ape::node.depth.edgelength(phy)[nodes],
              which.regimes = which.regimes))
}

weights <- function(a, lineage) {
  nt <- lineage$nodes_time
  ## Following is HUGE bottleneck, spends 60% ish time of a run with n=63
  # res <- c(vapply(head(seq_along(nt), -1),
  #                 function(i) exp(-a*(nt[1] - nt[i])) - exp(-a*(nt[1] - nt[i + 1])), FUN.VALUE = 0),
  #          exp(-a*nt[1]) ## Theta0\Ya
  #)
  res <- regimeWeightsHelper(a, nt) ## Rcpp wrapper, equivalent to above commented out code
  w <- vapply(lineage$which.regimes, function(e) sum(e*res), FUN.VALUE = 0)
  return(w)
}

weight.matrix <- function(phy, a, lineages){
  if(a > 300000000000) a <- 300000000000
  res <- t(vapply(lineages, function(x) weights(a, x), FUN.VALUE = numeric(length(lineages[[1]]$which.regimes))))
  #res <- t(sapply(lineages, function(x) weights(a, x)))
  rownames(res) <- phy$tip.label
  return(res)
}

## Thanks to user "snaut" at stackoverflow, http://stackoverflow.com/users/1999873/snaut
concat.factor <- function(...){
  as.factor(do.call(c, lapply(list(...), as.character)))
}