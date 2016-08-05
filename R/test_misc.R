## Accesories for tests in tests/testthat/
## Specify common parameters for all tests

#' @export
testsetup <- function(){
  pars <- list(myfunc = model.fit, 
               dummydata = read.table("dummy_data_slouch.txt"), 
               model.id.btk = "rReg")
  return(pars)
}



## MOdel types
## rReg
## IntcptReg