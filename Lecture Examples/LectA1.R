setwd("C:\\courses\\FISH 559_22\\TMB Workshop\\Lecture Examples\\")
data <- list(x=rivers)
parameters <- list(mu=0,logSigma=0)

require(TMB)
compile('LectA1A.cpp', flags="-Wno-ignored-attributes")
dyn.load(dynlib('LectA1A'))

##################
model <- MakeADFun(data,parameters,silent=T)
fit   <- nlminb(model$par, model$fn, model$gr)
rep   <- sdreport(model)
print(summary(rep))
model_tmp <- model
Function_num <- 0

model$fn_orig <- model$fn
model$fn <- function(x)
 {
  Function_num <<- Function_num + 1
  vv <- model$fn_orig(x)
  cat("fn",Function_num,vv,"\n")
  if (Function_num==30) {
    print(x)
    print(model$report(x))
    yy <- model$report(x)
   }
  #if (Function_num==30) AAA
  return(vv)
 }
model$gr_orig <- model$gr
model$gr <- function(x)
{
  vv <- model$gr_orig(x)
  cat("gr","Gradient",vv,"\n")
  if (Function_num==30) {
    print("Gradient")
    print(vv)
  }
  return(vv)
}
fit   <- nlminb(model$par, model$fn, model$gr)
rep   <- sdreport(model)
print(summary(rep))
