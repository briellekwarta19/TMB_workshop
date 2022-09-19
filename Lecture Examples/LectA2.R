
#setwd('Lecture Examples')
getwd()
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


#### maps ####
map <- list(u = c(factor(NA), factor(1), factor(NA), factor(2)))
#dont estimate 1st estimate 2nd, dont estimate 3rd, estimate 4th. 