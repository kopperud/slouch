## Semi-automatic test of dummy data
## Run with devtools::test() or ctrl + shift + T in Rstudio.

testpar <- testsetup()
list2env(testpar, envir = environment())
attach(dummydata)

## Check for skip
if (model.id.btk != "IntcptReg"){
  skip("Skipping IntcptReg")
}



## Testing intercept regression -- ultrametric
myfunc(ancestor, 
          time_ult,
          half_life_values = seq(0,45, length.out = 20),
          vy_values = seq(0,1.2, length.out = 35),
          response = trait_1,
          me.response = trait_1)

## Test intercept -- nonultrametric

myfunc(ancestor, 
          time,
          half_life_values = seq(0,45, length.out = 20),
          vy_values = seq(0,1.2, length.out = 35),
          response = trait_1,
          me.response = trait_1)
