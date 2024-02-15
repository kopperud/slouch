## Test environments
* Local OpenSuse, R 4.3.2
* Ubuntu Linux 20.04.1 LTS, R-release, GCC
* Fedora Linux, R-devel, clang, gfortran
* WinBuilder, R-devel (unstable) (2024-02-14 r85901 ucrt)

## R CMD check results

Appears go be all good

## Vignette caveats:
Vignettes may not build without R-packages "ape", "bookdown", and LaTeX packages "inputenc", "amsmath", "float". Vignettes were built without problems on winBuilder service, but a standard Windows install without LaTeX will have trouble.
