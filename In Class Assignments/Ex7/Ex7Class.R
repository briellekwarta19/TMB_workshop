library(TMB)

setwd("In Class Assignments\\Ex7")

Data <- read.csv("Ex7.csv")[,-1]
Nbatch <- max(Data$Batches)
Nunit <- length(Data[,1])

# Compile
require(TMB)
compile("Ex7Class.cpp", flags="-Wno-ignored-attributes")
dyn.load(dynlib("Ex7Class"))

# Create the data and parameters
data <- NULL
data$Nbatch=Nbatch; data$Nunit=Nunit;data$IsRandom=0
data$Batches <- Data$Batches; data$Units <- Data$Units; data$Treat <- Data$Treat 
data$Original <- Data$Original; data$Final <- Data$Final
#print(str(data))
parameters <- list(Control=0,Treatment=0,EpsB=rep(0,Nbatch),EpsU=rep(0,Nunit),logSigmaB=-1,logSigmaU=-1)
#print(str(parameters))


print("Model 1: No random effects")
map <- list(EpsB=rep(factor(NA),Nbatch),EpsU=rep(factor(NA),Nunit),logSigmaB=factor(NA),logSigmaU=factor(NA))

data$IsRandom <- 1

model <- MakeADFun(data, parameters,DLL="Ex7Class",map=map,control=list(eval.max=10000,iter.max=1000,rel.tol=1e-15),silent=T)
fit <- nlminb(model$par, model$fn, model$gr)
print(model$fn())
rep <- sdreport(model)
print(summary(rep))
Report1 <- model$report()
#print(Report)

print("Model 2: With random effects")
map <- NULL

data$IsRandom <- 2

model <- MakeADFun(data, parameters,DLL="Ex7Class", random = c("EpsB", "EpsU"),
                   map=map,control=list(eval.max=10000,iter.max=1000,rel.tol=1e-15),silent=T)
fit <- nlminb(model$par, model$fn, model$gr)
print(model$fn())
rep <- sdreport(model)
rep2 <- summary(rep)
print(rep2[row.names(rep2)!="EpsU" & row.names(rep2)!="EpsB",])
Report2 <- model$report()
#print(Report)

#Pearson residuals
data_p <- data$Final/data$Original

Fixed_resids <- (data_p - Report1$Prob)/ (Report1$Prob*(1-Report1$Prob))
Random_resids <- (data_p - Report2$Prob)/ (Report2$Prob*(1-Report2$Prob))

par(mfrow = c(2,1))
hist(Fixed_resids)
hist(Random_resids)
