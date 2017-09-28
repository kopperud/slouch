## Tests weightmatrix
concat.factor <- function(...){
  as.factor(do.call(c, lapply(list(...), as.character)))
}

library(ape)
set.seed(15)
n <- 85

phy <- rtree(n)
regimes_tip <- sample(c("A", "B"), n, replace = TRUE)

ans <- ace(regimes_tip, phy, type = "discrete")
regimes_internal <- factor(colnames(ans$lik.anc))[apply(ans$lik.anc, 1, which.max)]

regimes <- slouch::concat.factor(regimes_tip, regimes_internal)

lineages <- lapply(1:n, function(e) lineage.constructor(phy, e, regimes))


alphas <- c(0, 1, 15, 30000000000)

naive_regweightshelper <- function(a, nt){
  res <- c(vapply(head(seq_along(nt), -1),
                  function(i) exp(-a*(nt[1] - nt[i])) - exp(-a*(nt[1] - nt[i + 1])), FUN.VALUE = 0),
           exp(-a*nt[1]) ## Theta0\Ya
  )
  return(res)
}

#slouch::weightshelper(3, lineages[[35]]$nodes_time)

#regimeweights <- naive_regweightshelper(0, nt)

for (a in alphas){
  old <- lapply(lineages, function(lineage) naive_regweightshelper(a, lineage$nodes_time))
  new <- lapply(lineages, function(lineage) weights_segments(a, lineage))
  # print(microbenchmark::microbenchmark(old <- lapply(lineages, function(lineage) naive_regweightshelper(a, lineage$nodes_time)),
  #                                      new <- lapply(lineages, function(lineage) weights_segments(a, lineage)),
  #                                      rcpp <- lapply(lineages, function(lineage) regimeWeightsHelper(a, lineage$nodes_time))))
  # 
  test_that(paste("Regimeweights for alpha =", a), {
    expect_equal(old, new)
  })
  test_that(paste("Coef sum to one for alpha =", a), {
    sum_coef <- lapply(new, sum)
    expect_equal(sum_coef, as.list(rep(1,n)))
  })
}

