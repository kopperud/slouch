---
output: html_document
---
# How to: Install & load devtools

Careful, it is in development.
```
install.packages("devtools")
library(devtools)
```

# Install and load this packaage
```
devtools::install_github("bjornkopperud/slouchexp")
library(slouchexp)
```

Since it is in development, I've changed the name of model.fit command, as to point out that it is in dev, to `model.fit.dev2()`.

Rudimentary helpfiles exist, such as `?model.fit.dev2`
