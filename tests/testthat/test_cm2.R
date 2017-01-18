## Formal test to see whether vectorized implementation of cm2, part of residual 
## var-cov matrix is computationally equivalent to the old, slow loop version

# ## function to generate arbitrary posdef matrix
# foo <- function(m){
#   crossprod(matrix(abs(rnorm(m*m)), ncol=m))
# }

library(ape)
N <- ceiling(runif(1, 30, 50))
phy <- rtree(N)
a <- 3

mrca1 <- mrca(phy)
times <- node.depth.edgelength(phy)
ta <- matrix(times[mrca1], nrow=N, dimnames = list(phy$tip.label, phy$tip.label))
T.term <- times[1:N]
tia <- times[1:N] - ta
tja <- t(tia)
tij <- tja + tia

## Old calculation
cm2 <- matrix(0, ncol=N, nrow=N)
num.prob <- matrix(0, ncol=N, nrow=N)

for(p in 1:N)
{
  for(q in 1:N)
  {
    if(ta[p,q]==0)num.prob[p,q]=1 else num.prob[p,q]=(1-exp(-a*ta[p,q]))/(a*ta[p,q]);
  }
}
for(p in 1:N)
{
  for(q in 1:N)
  {
    cm2[p,q]<-(((1-exp(-a*T.term[p]))/(a*T.term[p]))*((1-exp(-a*T.term[q]))/(a*T.term[q]))-(exp(-a*tia[p, q])*(1-exp(-a*T.term[q]))/ (a*T.term[q])+ exp(-a*tja[p, q])*(1-exp(-a*T.term[p]))/(a*T.term[p]))*(num.prob[p,q]));
  }
}

  
## New version
cm2_vectorized = unname(calc.cm2(a, T.term, N, tia, tja, ta))

test_that("residual var-cov matrix is identical in vectorized form", {
  expect_equal(cm2, cm2_vectorized)
})