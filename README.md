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
