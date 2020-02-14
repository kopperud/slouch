## Test environments
* Local Debian 10, R 3.5.22
* Ubuntu 16.04.6 LTS (on Travis-ci), R 3.6.2
* WinBuilder R-devel (unstable) (2020-01-28 r77738)

## R CMD check results
Status: 1 NOTE

Changed maintainer e-mail in DESCRIPTION. The words "Ornstein" and "Uhlenbeck" are surnames, and are not misspelled.


## Vignette caveats:
Vignettes may not build without R-packages "ape", "bookdown", and LaTeX packages "inputenc", "amsmath", "float". Vignettes were built without problems on winBuilder service, but a standard Windows install without LaTeX will have trouble.
