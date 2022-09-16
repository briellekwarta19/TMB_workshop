setwd("C:\\courses\\FISH 559_22\\TMB Workshop\\Lecture Examples\\")
data <- read.table("LectA2.dat", header=TRUE)
parameters <- list(b0=0, b1=0, logSigma=0)

require(TMB)
compile("LectA2.cpp", flags="-Wno-ignored-attributes")
dyn.load(dynlib("LectA2"))

################################################################################

model <- MakeADFun(data, parameters, DLL="LectA2",silent=T)
fit <- nlminb(model$par, model$fn, model$gr)

best <- model$env$last.par.best
rep <- sdreport(model)

print(best)
print(rep)

print(data$y)
simdata <- model$simulate(complete = TRUE)
print(simdata$y)
simdata <- model$simulate(complete = TRUE)
print(simdata$y)