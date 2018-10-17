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


