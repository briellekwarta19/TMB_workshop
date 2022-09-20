setwd("In Class Assignments/Ex2")

data <- read.table("EX2.dat", header=TRUE)

require(TMB)
# compile("Example2.cpp", flags="-Wno-ignored-attributes")
# dyn.load(dynlib("Example2"))

compile("Ex2.cpp", flags="-Wno-ignored-attributes")
dyn.load(dynlib("Ex2"))

################################################################################
#### Model 1 ####
NP <- as.matrix(data[,2:4])
rownames(NP) <- NULL
PY <- as.matrix(data[,5:7])
rownames(PY) <- NULL

data_mod1 <- list(NC = as.vector(data$Predator), NP = NP, PY= PY , Model = 1)
map1 <- list(a0 = factor(NA), logk = factor(NA))

parameter1 <- list(loga=rep(1,3))

model <- MakeADFun(data_mod1, 
                   parameter1,
                   DLL="Ex2")

print(attributes(model))

fit <- nlminb(model$par, model$fn, model$gr)
for (i in 1:3)
 fit <- nlminb(model$env$last.par.best, model$fn, model$gr)
              #STOP AT THE BEST VALUE, then start again at best value

rep <- sdreport(model)

# Sumamrize ALL
print(summary(rep,p.value=T))
# Restrict what comes out to the fixed parameters only
print(summary(rep,select="fixed",p.value=F))

#### Model 2 ####