[![Build Status](https://travis-ci.org/kopperud/slouch.svg?branch=master)](https://travis-ci.org/kopperud/slouch) [![codecov.io](https://codecov.io/github/kopperud/slouch/coverage.svg?branch=master)](https://codecov.io/github/kopperud/slouch?branch=master)

# SLOUCH: Stochastic Linear Ornstein-Uhlenbeck Comparative Hypotheses

SLOUCH is an R-package implementation of some statistical models that are used in evolutionary biology, more specifically phylogenetic comparative methods. Given a phylogenetic tree and a comparative dataset, SLOUCH allows the user to fit models that are conceptually consistent with adaptation towards an optimal state, and includes options to incorporate observational or measurement error in regression analyses. See Hansen *et al.* (2008) for the original presentation of SLOUCH.

# Install instructions

The R-package `devtools` makes it easy to install R packages straight from github.
```
install.packages("devtools")
library(devtools)

devtools::install_github("kopperud/slouch")
library(slouch)
```


# Documentation

Standard R-package documentation can be seen by entering the command `?slouch`, or visiting the package website [the package website](https://kopperud.github.io/slouch/).


### References

* Hansen, T. F., Pienaar, J., & Orzack, S. H. (2008). A comparative method for studying adaptation to a randomly evolving environment. Evolution, 62(8), 1965-1977.