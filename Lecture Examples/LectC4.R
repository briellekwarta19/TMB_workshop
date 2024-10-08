setwd("Lecture Examples")
data <- read.table("LectA2.dat", header=TRUE)
parameters <- list(b0=0, b1=0, logSigma=0)

require(TMB)
compile("LectA2.cpp", flags="-Wno-ignored-attributes")
dyn.load(dynlib("LectA2"))

################################################################################


print("Phase 1; No log-Sigma")
map <- list(logSigma=factor(NA))
model <- MakeADFun(data, parameters, DLL="LectA2",map=map,silent=T)
fit <- nlminb(model$par, model$fn, model$gr)
best <- model$env$last.par.best
print(as.numeric(best))

print("Phase 2; All parameters - note I ordered my parameters to make this easier")
Phase2Init <- c(as.numeric(best),0)
model <- MakeADFun(data, parameters, DLL="LectA2",silent=T)
fit <- nlminb(Phase2Init, model$fn, model$gr)
best <- model$env$last.par.best
rep <- sdreport(model)
print(best)
print(rep)

