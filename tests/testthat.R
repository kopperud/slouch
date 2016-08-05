library(testthat)
library(slouch)

myfunc <- model.fit.dev ## use model.fit to test old source.R || model.fit.dev to test package
#dummydata <- read.table("/home/bjorn/Downloads/dummy_data_slouch.txt")
dummydata <- read.table("dummy_data_slouch.txt")

## Test only one function at a time
model.id.btk <- "IntcptReg"


testsetup <- function(){
  pars <- list(myfunc = myfunc, dummydata = dummydata, model.id.btk = model.id.btk)
  return(pars)
}

test_check("slouch")

