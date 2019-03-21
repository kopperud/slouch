context("ANOVA (diet)")
## Generate tree for testing test
library(ape)


## Simulation, test broken due to R 3.6.0 release changing the sample 
## method, and thus `set.seed(4)` is broken.
  
# n = 80
# set.seed(4)
# phy <- rtree(n)
# trait_1 <- runif(n)
# trait_1_SE_sq <- 0.01*runif(n)
# regimes_tip <- factor(sample(c("A", "B"), n, replace = TRUE))
# 
# ## Ancestral state recon
# ans <- ace(regimes_tip, phy, type ="d")
# regimes_internal <- factor(levels(regimes_tip)[apply(ans$lik.anc, 1, function(e) which.max(e))])
# phy$node.label <- regimes_internal
# 
# regimes <- concat.factor(regimes_tip, regimes_internal)
# lineages <- lapply(1:n, function(e) lineage.constructor(phy, e, regimes)) #; names(lineages) <- phy$tip.label
#
## Plot
#plot(phy); tiplabels(regimes_tip); nodelabels(regimes_internal)



data("neocortex")
data("artiodactyla")

neocortex <- neocortex[match(artiodactyla$tip.label, neocortex$species), ]

m1 <- slouch.fit(artiodactyla,
                     species = neocortex$species,
                     hl_values = seq(0,0.4, length.out = 15),
                     vy_values = seq(0.05,0.15, length.out = 15),
                     response = neocortex$neocortex_area_mm2_log_mean,
                     mv.response = neocortex$neocortex_se_squared,
                     fixed.fact = neocortex$diet)



reference <- structure(c(-28323.9802125233, 53231.086989994, 9.91679560722592, 
                         19774.1791676592, 30572.5096605727, 0.377428304377147), .Dim = 3:2, .Dimnames = list(
                           c("Br", "Gr", "MF"), c("Estimates", "Std. error")))


test_that("No errors ANOVA (neocortex diet)", {
  expect_equal(reference,
               m1$beta_primary$coefficients,
               tolerance = 9e-04)
})
