setwd("Lecture Examples")

library(TMB)
compile("LectC5.cpp",libinit=FALSE, flags="-Wno-ignored-attributes")                    ## notice flag
dyn.load(dynlib("LectC5"))
set.seed(123)
data <- list(Y = rnorm(10) + 1:10, x=1:10)
parameters <- list(a=0, b=0, logSigma=0)

#obj <- MakeADFun(data, parameters, DLL="LectD42")


## setup call from R
myline<-function(x,a,b).Call("call_myline",as.double(x), as.double(a), as.double(b),PACKAGE="LectC5")

## now you can call

print(myline(1:10, 1, 1.0001))

