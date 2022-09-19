library(here)
data <- read.table(here::here("In Class Assignments", "Ex1", "ex1.dat"), header=TRUE)
colnames(data) <- c("Age", "Length")
parameters <- list(a0=1, loga50=1, logSigma=1, logLinf = 1, logk = 1)

require(TMB)
compile("ex1TMB.cpp", flags="-Wno-ignored-attributes")
dyn.load(dynlib("ex1TMB"))

################################################################################
moddata_1 <- list(Age = data$Age, Length = data$Length, Model = 1)
map_1 <- list(a0 = factor(NA), logk = factor(NA))
model_1 <-MakeADFun(moddata_1, map = map_1, parameters, DLL="ex1TMB",silent=T)
fit_1 <- nlminb(model_1$par, model_1$fn, model_1$gr)

best_1 <- model_1$env$last.par.best
rep_1 <- sdreport(model_1)

print(best_1)
print(rep_1)


dim(data)

################################################################################
moddata_2 <- list(Age = data$Age, Length = data$Length, Model = 2)
map_2 <- list(loga50 = factor(NA))
model_2 <-MakeADFun(moddata_2, map = map_2, parameters, DLL="ex1TMB",silent=T)
fit_2 <- nlminb(model_2$par, model_2$fn, model_2$gr)

best_2 <- model_2$env$last.par.best
rep_2 <- sdreport(model_2)

print(best_2)
print(rep_2)

