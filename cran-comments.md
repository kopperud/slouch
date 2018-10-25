## Test environments
* Local Fedora 27, R 3.5.0
* Ubuntu 14.04.5 LTS (on Travis-ci), R 3.5.1

## R CMD check results
Status: OK

Apologies for fast resubmission, one vignette should be retracted. The words "Ornstein" and "Uhlenbeck" are surnames, and are not misspelled.


## Vignette caveats:
Vignettes may not build without R-packages "ape", "bookdown", and LaTeX packages "inputenc", "amsmath", "float". Vignettes were built without problems on winBuilder service, but a standard Windows install without LaTeX will have trouble.
