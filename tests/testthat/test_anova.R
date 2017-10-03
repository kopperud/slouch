
## Generate tree for testing test
library(ape)
n = 80
set.seed(4)
phy <- rtree(n)
trait_1 <- runif(n)
trait_1_SE_sq <- 0.01*runif(n)
#regimes_tip <- factor(sample(c("A", "B", "C"), n, replace = TRUE))
regimes_tip <- factor(sample(c("A", "B"), n, replace = TRUE))

## Ancestral state recon
ans <- ace(regimes_tip, phy, type ="d")
regimes_internal <- factor(levels(regimes_tip)[apply(ans$lik.anc, 1, function(e) which.max(e))])
phy$node.label <- regimes_internal
#levels(regimes_internal) <- c(levels(regimes_internal), "root"); regimes_internal[1] <- "root"

## Plot
#plot(phy); tiplabels(regimes_tip); nodelabels(regimes_internal)

regimes <- concat.factor(regimes_tip, regimes_internal)
lineages <- lapply(1:n, function(e) lineage.constructor(phy, e, regimes)) #; names(lineages) <- phy$tip.label


m1 <- slouch.fit(phy,
                     species = phy$tip.label,
                     hl_values = seq(0,0.4, length.out = 15),
                     vy_values = seq(0.05,0.15, length.out = 15),
                     response = trait_1,
                     me.response = trait_1_SE_sq,
                     fixed.fact = regimes_tip)


#detach(dummydata)

test_that("No errors ANOVA", {
  expect_equal(structure(c(0.421599539503631, 0.481516052311344, 0.0453131399491333, 
                           0.050412202001475), .Dim = c(2L, 2L), .Dimnames = list(c("A", 
                                                                                    "B"), c("Estimates", "Std. error"))),
               m1$opt.reg$coefficients,
               tolerance = 9e-07)
})
