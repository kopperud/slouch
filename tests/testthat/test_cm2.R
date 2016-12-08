## Formal test to see whether vectorized implementation of cm2, part of residual 
## var-cov matrix is computationally equivalent to the old, slow loop version

## function to generate arbitrary\force symmetry on matrix
foo <- function(m){
  m[lower.tri(m)] = t(m)[lower.tri(m)]
  return(m)
}
N <- 30
a <- 0.56
tia <- foo(matrix(abs(rnorm(N*N)), ncol=N))
tja <- foo(matrix(abs(rnorm(N*N)), ncol=N))
ta <- foo(matrix(abs(rnorm(N*N)), ncol=N))
#T.term <- matrix(rep(1, N), ncol=1)
T.term <- matrix(rnorm(N), ncol=1)

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
cm2_vectorized = calc.cm2(a, T.term, N, tia, tja, ta)

test_that("residual var-cov matrix is identical in vectorized form", {
  expect_equal(cm2, cm2_vectorized)
})