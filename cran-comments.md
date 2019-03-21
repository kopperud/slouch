## Test environments
* Local Fedora 27, R 3.5.22
* Ubuntu 14.04.5 LTS (on Travis-ci), R 3.5.2
* WinBuilder R 3.6.0 (Devel)

## R CMD check results
Status: OK

Amended tests not to depend on `set.seed()`, such that it works with new sampler in R 3.6.0. The words "Ornstein" and "Uhlenbeck" are surnames, and are not misspelled.


## Vignette caveats:
Vignettes may not build without R-packages "ape", "bookdown", and LaTeX packages "inputenc", "amsmath", "float". Vignettes were built without problems on winBuilder service, but a standard Windows install without LaTeX will have trouble.
