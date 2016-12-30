## Test run time dummy dataset


dummydata <- read.table("dummy_data_slouch.txt")


## Square the SE, so that it is the Var of the mean
dummydata$SE_sq_trait_1 <- NA
dummydata$SE_sq_trait_1[!is.na(dummydata$SE_trait_1)] <- (dummydata$SE_trait_1[!is.na(dummydata$SE_trait_1)])^2

dummydata$SE_sq_trait_2 <- NA
dummydata$SE_sq_trait_2[!is.na(dummydata$SE_trait_2)] <- (dummydata$SE_trait_1[!is.na(dummydata$SE_trait_2)])^2

## Fake data, rnorm
dummydata$fake_rnorm <- NA
dummydata$fake_rnorm[!is.na(dummydata$species)] <- rnorm(length(na.exclude(dummydata$species)))

attach(dummydata)

m1 <- model.fit.dev2(dummydata$ancestor,
                     dummydata$time,
                     half_life_values = seq(0.1,5, length.out = 3),
                     vy_values = seq(0.4,7, length.out = 3),
                     response = dummydata$trait_1,
                     me.response = dummydata$SE_sq_trait_1,
                     fixed.fact = dummydata$three_niches,
                     ultrametric = TRUE)

detach(dummydata)

test_that("No errors ANOVA", {
  expect_equal(structure(c(1.72188320563384, 0.667703233610768, 0.951808167895859, 
                           1.2482851617156, 1.06450145357172, 0.892032143040452), .Dim = c(3L, 
                                                                                           2L), .Dimnames = list(c("ANC", "A", "B"), c("Estimates", "Std. error"
                                                                                           ))),
               m1$opt.reg$coefficients)
})