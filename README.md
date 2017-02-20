---
output:
  html_document: default
  pdf_document: default
  word_document: default
---
[![Build Status](https://travis-ci.org/kopperud/slouch.svg?branch=master)](https://travis-ci.org/kopperud/slouch)

# Install and load devtools

Devtools makes it much easier to install R packages straight from github
```
install.packages("devtools")
library(devtools)
```

# Install and load slouch
```
devtools::install_github("kopperud/slouch")
library(slouch)
```

Since it is in development, I've changed the name of model.fit command, as to point out that it is in dev, to `model.fit.dev2()`.

Rudimentary helpfiles exist, such as `?model.fit.dev2`

# Notable changes since v2.*

### Bugfixes:
* Regression coefficients now correctly start with an OLS estimate instead of using the GLS-estimate for the previous iteration in the gridsearch, for categorical models
* For categorical models: Seeding OLS regression is now correctly done for every iteration in grid search, for every model matrix X, instead of an invariate model matrix X whose alpha is arbitrary (e.g a = 10)
* Most if not all model combinations should now support multiple continuous covariates, e.g fixed.cov = cbind(X, Z, ...), without crashing before returning model outputs.
* Vectorize the code calculating residual var-cov matrix for models with brownian covariates, thus fixing one index typo. Bug was only relevant for non-ultrametric trees.
* Fix incorrect calculation of branch lengths "tja" and "tia" for non-ultrametric trees, see below for new tree format
* The weight matrix generated with categorical covariates is no longer sorted by columns, the order is now invariant with respect to alpha. This was not strictly a bug but if the coefficient labels are not rearranged in the same way, output will be wrong & misleading.

### Structural changes
* Code duplication largely removed, R-code reduced from about 4k to 1k lines.
* SLOUCH now uses Git version control.
* Direct effect models now about 4x faster. Models with categorical and/or brownian covariates are now estimated nearly as fast as direct effect models.
* Grid-search now has support for using multiple CPU cores, currently only for windows.
* A hillclimber routine (L-BFGS) implemented in R base `optim()` is implemented. This routine can be used standalone or combined with the hillclimber, where the hillclimber starts at the best ML-estimate from the grid-search.
* The format used to display trees in OUCH/SLOUCH is now scrapped in favor of package `ape`. Calculation of branch lengths is now dependent on this package, and the weight matrix functions have been rewritten from scratch to accomodate the new format. SLOUCH is now technically no longer an extension of OUCH.
* Inner "workhorse" loop of the weight matrix functions is rewritten in C++
* Because of this, functions for simulating OU traits, plotting the tree, and the Fitch algorithm for reconstructing ancestral states have been scrapped.
* SLOUCH now depends on packages `parallel`, `ape`, `Rcpp`.
* Model outputs are no longer printed directly in the console, but returned as a composite object of class `slouch`, essentially a list with lots of information and methods for printing and plotting.
* Arguments `intercept` and `ultrametric` are misleading and to be changed.
* SLOUCH now has one unit test, should be many more

### Other changes
* For grid-search, the marginal support set for alpha and sigma are now calculated and given in a table.
* The GLS estimator is now not used as-is to calculate regression coefficients. Instead, the lower triangular matrix of the cholesky factorization of $V = C^TC$, $L = C^T$ is used to "transform" $X_{*} = L^{-1}X$ and $Y_{*} = L^{-1}Y$. Next, $X_*$ and $Y_*$ are fed to R's `lm.fit()`, which uses QR-decomposition to calculate coefficients. This makes for about a 4x speedup. Also if X is singular, the program will now not crash but return "NA" in the regression coefficients.
* Calculation of $log(det(V))$ is now done differently, both to increase speed and numerical stability (not returning `Inf`) when N is large.
* Docstrings are now available, though crude.

