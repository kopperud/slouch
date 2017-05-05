## ---- fig.width=7--------------------------------------------------------
library(ape)

set.seed(8)
n <- 80
phy <- rtree(n)

oldmar <- par("mar"); par(mar = c(0,0,0,0))
plot(ladderize(phy))
par(mar = oldmar)

## ---- fig.width = 7------------------------------------------------------
mydata <- data.frame(species = phy$tip.label,
                     y = rnorm(n, sd = 4) + 1:n*1.7,
                     x = rnorm(n, sd = 4) + 1:n,
                     z = rnorm(n, sd = 4) + 1:n*(-0.6))

plot(mydata$x, mydata$y, ylab = "y", xlab = "x")

## Scramble the rows
mydata <- mydata[sample(1:n, n),]

## Check whether they are lined up correctly
mydata$species == phy$tip.label

## ------------------------------------------------------------------------
mydata <- mydata[match(phy$tip.label, mydata$species) ,]

## Check if they line up again
mydata$species == phy$tip.label

## ---- fig.show='hold', fig.width=3, fig.height=3-------------------------
library(slouch)
model0 <- slouch.fit(phy = phy,
                     hl_values = seq(0, 50, length.out = 10),
                     vy_values = seq(150, 800, length.out = 10),
                     species = mydata$species,
                     response = mydata$y)
plot(model0, cex.lab = 0.7, cex.axis = 0.7)

model0 <- slouch.fit(phy = phy,
                     hl_values = seq(0, 100, length.out = 25),
                     vy_values = seq(0.1, 10000, length.out = 25),
                     species = mydata$species,
                     response = mydata$y)
plot(model0, cex.lab = 0.7, cex.axis = 0.7)

## ------------------------------------------------------------------------
model0 <- slouch.fit(phy = phy,
                     species = mydata$species,
                     response = mydata$y,
                     hillclimb = TRUE)
plot(model0)
model0$oupar

## ---- fig.show='hold'----------------------------------------------------
model1 <- slouch.fit(phy = phy,
                     species = mydata$species,
                     response = mydata$y,
                     fixed.cov = mydata$x,
                     lower = c(0, 0.01),
                     hillclimb = TRUE)
plot(model1)

plot(mydata$x, mydata$y, ylab = "y", xlab = "x", main = "Trait plot")
abline(lm(mydata$y ~ mydata$x), col = "black", lwd = 2)
abline(model1$opt.reg$coefficients[,1], col = "orange", lwd = 2)
model1$opt.reg$coefficients

## ------------------------------------------------------------------------
model2 <- slouch.fit(phy = phy,
                     species = mydata$species,
                     response = mydata$y,
                     fixed.cov = cbind(x = mydata$x, z = mydata$z),
                     lower = c(0, 0.01),
                     hillclimb = TRUE)
plot(model2)

## ------------------------------------------------------------------------
categories <- c("A", "B", "C")
mydata$category <- sample(categories, n, replace = TRUE)

## ---- fig.width=7--------------------------------------------------------
library(ape)
reconstruction <- ace(mydata$category, phy, type = "d")

## Extract the most likely regime for each internal node
## These have order n+1, n+2, n+3 ...
internal_regimes <- apply(reconstruction$lik.anc, 
                          1, 
                          function(e) colnames(reconstruction$lik.anc)[which.max(e)])

## Concatenate tip and internal regimes. These will have order 1,2,3,4, ...
regimes <- c(mydata$category, internal_regimes)

## Pick out the regimes of the edges, in the order of phy$edge
edge_regimes <- factor(regimes[phy$edge[,2]])

oldmar <- par("mar"); par(mar = c(0,0,0,0))
plot(phy, edge.color = c("Black", "#EE7600", "blue")[edge_regimes], edge.width = 3)
par <- par(oldmar)

## ------------------------------------------------------------------------
phy$node.label <- internal_regimes

model3 <- slouch.fit(phy = phy,
                     species = mydata$species,
                     response = mydata$y,
                     fixed.cov = mydata$x,
                     fixed.fact = mydata$category,
                     hillclimb = TRUE,
                     lower = c(0, 0.01))

model3$opt.reg$coefficients

## ---- fig.width = 7------------------------------------------------------
phy$node.label <- internal_regimes

model4 <- slouch.fit(phy = phy,
                     hl_values = seq(0, 0.25, length.out = 25),
                     vy_values = seq(20, 60, length.out = 25),
                     species = mydata$species,
                     response = mydata$y,
                     fixed.cov = mydata$x,
                     random.cov = mydata$z,
                     fixed.fact = mydata$category,
                     hillclimb = TRUE,
                     lower = c(0, 0.1))

plot(model4)
model4$opt.reg$coefficients

## ---- fig.width = 7------------------------------------------------------
phy$node.label <- internal_regimes

model5 <- slouch.fit(phy = phy,
                     hl_values = seq(0, 0.25, length.out = 25),
                     vy_values = seq(20, 60, length.out = 25),
                     species = mydata$species,
                     response = mydata$y,
                     me.response = 0.01*rnorm(n),
                     fixed.cov = mydata$x,
                     me.fixed.cov = 0.01*rnorm(n),
                     random.cov = mydata$z,
                     me.random.cov = 0.01*rnorm(n),
                     hillclimb = TRUE,
                     lower = c(0, 0.1))

plot(model5)

