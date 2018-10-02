## Test environments
* Local Fedora 27, R 3.5.0
* Ubuntu 14.04.5 LTS (on Travis-ci), R 3.5.1
* winBuilder x86_64-w64-mingw32, R 3.5.1

## R CMD check results
0 errors | 0 warnings | 1 notes

The note says that a URL is invalid, and possibly some mis-spelled words. I can visit the URL just fine (after redirect), don't know why there is an error. The words "Ornstein" and "Uhlenbeck" are surnames, and are not misspelled.

I read the CRAN URL checks
https://cran.r-project.org/web/packages/URL_checks.html


```
	Possibly mis-spelled words in DESCRIPTION:
	  Ornstein (3:26, 12:106)
	  Uhlenbeck (3:35, 12:115)

	Found the following (possibly) invalid URLs:
	  URL: https://doi.org/10.2307/2411186
	    From: inst/doc/background.html
	          inst/doc/examples.html
	    Status: 403
	    Message: Forbidden
```

## Vignette caveats:
Vignettes may not build without R-packages "ape", "bookdown", and LaTeX packages "inputenc", "amsmath", "float". Vignettes are built without problems on winBuilder service, but a standard Windows install without LaTeX will have trouble.
