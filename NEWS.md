# slouch 2.1.5

* Include functionality for allowing SIMMAP style trees to dictate the regimes painted on the tree. In other words, also allowing the regimes to change character state along the branches, instead of just at the branching events. 
* Fixed bug where if there were multiple predictor variables modeled as BM, then the predictors were re-ordered
* Bugfix for allowing a Brownian-motion model with a single trend (i.e. estimating the ancestral state and a single trend parameter, but no regime map)
* Don't print the trend contrasts in the summary unless those transitions were actually present on the tree

# slouch 2.1.4

* Changed maintainer e-mail
* Bugfix: model with random, continuous covariate was broken. Bug only relevant for installs from github, between commits ee5d38 and 04b103. CRAN repository not affected.

# slouch 2.1.3

* Added rudimentary support for interaction terms between direct-effect continuous predictors, and categorical predictors. Caveat: measurement error in continuous predictor not supported for estimating interaction terms.

# slouch 2.1.2

* Amended tests not to depend on `set.seed()`, such that it works with new sampler in R 3.6.0. 

# slouch 2.1.1

* Removed 'background' vignette

# slouch 2.1.0

* Added a `NEWS.md` file to track changes to the package.
* Breaking changes: several variable names in function arguments and model outputs have been changed for clarity.
* Documentation improved.
* Model summary print changed
* Draft for long-form manual/vignette begun.

# slouch 2.0.0

slouch 2.0.0 marks a near-complete rewrite of the entire package.

## New features

* Slouch now uses the phylogenetic tree format from package `ape`.
* Model outputs are returned as a composite object, which consequentely is programmable.
* Implementation of explicit Brownian-motion models, including intercept-only models, regression with direct-effect covariates, random-effect (trend) covariates, regime-dependent trends and options to estimate `Ya`.
* Slouch now allows estimation of parameters using numerical optimization techniques.
* Computation in general is much faster, due to better memory management and use of Cholesky decomposition in regression coefficient estimation.
* Docstrings begun.

## Bug fixes

* Fixed several bugs when calculating the residual variance-covariance matrix, in particular for random-effect models with non-ultrametric trees.
* Fixed incorrect calculation of `tja` and `tia` for non-ultrametric trees.

## Scrapped features

* Fitch algorithm for ancestral state character estimation.
* The `slouchtree` or `ouchtree` phylogeny format.


