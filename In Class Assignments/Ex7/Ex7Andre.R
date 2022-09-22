library(TMB)

setwd("C:\\courses\\FISH 559_22\\TMB Workshop\\In Class Assignments\\Ex7")

Data <- read.csv("Ex7.csv")[,-1]
Nbatch <- max(Data$Batches)
Nunit <- length(Data[,1])

# Compile
require(TMB)
compile("Ex7.cpp", flags="-Wno-ignored-attributes")
dyn.load(dynlib("Ex7"))

# Create the data and parameters
data <- NULL
data$Nbatch=Nbatch; data$Nunit=Nunit;data$IsRandom=0
data$Batches <- Data$Batches; data$Units <- Data$Units; data$Treat <- Data$Treat 
data$Original <- Data$Original; data$Final <- Data$Final
#print(str(data))
parameters <- list(Control=0,Treatment=0,EpsB=rep(0,Nbatch),EpsU=rep(0,Nunit),logSigmaB=-1,logSigmaU=-1)
#print(str(parameters))

print("Model 1: No random effects")
data$IsRandom <- 0
map <- list(EpsB=rep(factor(NA),Nbatch),EpsU=rep(factor(NA),Nunit),logSigmaB=factor(NA),logSigmaU=factor(NA))
model <- MakeADFun(data, parameters,DLL="Ex7",map=map,control=list(eval.max=10000,iter.max=1000,rel.tol=1e-15),silent=T)
fit <- nlminb(model$par, model$fn, model$gr)
print(model$fn())
rep <- sdreport(model)
print(summary(rep))
Report1 <- model$report()
cat("Objective fn",Report1$f,"Negative logL",Report1$g,"; deviance:",Report1$Deviance,"\n")
#print(Report)

print("Model 2: With random effects")
#map <- list(logSigmaB=factor(NA),logSigmaU=factor(NA),EpsU=rep(factor(NA),Nunit))
map <- NULL
data$IsRandom <- 1
model <- MakeADFun(data, parameters,DLL="Ex7",map=map,random=c("EpsU","EpsB"),control=list(eval.max=10000,iter.max=1000,rel.tol=1e-15),silent=T)
fit <- nlminb(model$par, model$fn, model$gr)
print(model$fn())
rep <- sdreport(model)
rep2 <- summary(rep)
print(rep2[row.names(rep2)!="EpsU" & row.names(rep2)!="EpsB",])
Report2 <- model$report()
cat("Objective fn",Report2$f,"Negative logL",Report2$g,"; deviance:",Report2$Deviance,"\n")
#print(Report)

# Residuals
par(mfrow=c(2,2))
Resids <- (data$Final/data$Original-Report1$Prob)/(Report2$Prob*(1.0-Report2$Prob))
hist(Resids)
plot(1:length(Resids),Resids,xlab="Data point")
Resids <- (data$Final/data$Original-Report2$Prob)/(Report2$Prob*(1.0-Report2$Prob))
hist(Resids)
plot(1:length(Resids),Resids,xlab="Data point")

