
################################################################################

setwd("In Class Assignments\\Ex5")
require(TMB)

################################################################################

m <- scan("ex5.dat",skip=1,n=1,quiet=T)
TT <- scan("ex5.dat",skip=3,n=m,quiet=T)
Tmax <- scan("ex5.dat",skip=5,n=1,quiet=T)
B <- matrix(scan("ex5.dat",skip=7,n=m*Tmax,quiet=T),ncol=Tmax,byrow=T)
R <- matrix(scan("ex5.dat",skip=58,n=m*Tmax,quiet=T),ncol=Tmax,byrow=T)
Phi0 <- scan("ex5.dat",skip=109,n=m,quiet=T)

data <- list(m=m,TT=TT,Tmax=Tmax,B=B,R=R,Phi0=Phi0)

################################################################################
# Hint Nproj = 0 for basic estimation
# Basic estimation
################################################################################

# Set the parameters to initial values
#bounds are on ex 5 document
# parameters <- list(dummy=0,mu=-1,log_tau=2,B0=rep(500, 50),
#                    log_sigR=rep(log(2),50),eta=rep(0.2,50))


parameters <- list(dummy=0,mu=0,log_tau=-1.5,B0=rep(1000.1/2, 50),
                   log_sigR=rep(-1.5,50),eta=rep(0,50))


#print(data)
#print(parameters)

compile("Ex5Class.cpp", flags="-Wno-ignored-attributes")
dyn.load(dynlib("Ex5Class"))

# test code - for checking for minimization
map<-list(mu=factor(NA),log_tau=factor(NA),B0=rep(factor(NA),m),log_sigR=rep(factor(NA),m),eta=rep(factor(NA),m))
model <- MakeADFun(data, parameters, DLL="Ex5Class",silent=T,map=map)
xx <- model$fn(model$env$last.par)
cat(xx,model$report()$obj_fun,"\n")
#AAA

# Now estimate everything
map<-list(dummy=factor(NA))
model <- MakeADFun(data, parameters, random="eta", DLL="Ex5Class",silent=T,map=map)

rep <- sdreport(model)
rep
# Bounds on the parameters
# parameters <- list(dummy=0,mu=-1,log_tau=2,B0=rep(500, 50),
#                    log_sigR=rep(log(2),50),eta=rep(5,50))

lowbnd= c(mu= -10, log_tau = -7, B0= rep(0.1,50),
         log_sigR = rep(-7, 50), eta= rep(-Inf, 50))

uppbnd= c(mu =  10, log_tau = 4, B0 = rep(1000,50),
          log_sigR = rep(4, 50), eta = rep(Inf, 50))

fit <- nlminb(model$par, model$fn, model$gr,
              control=list(rel.tol=1e-12,eval.max=100000,iter.max=1000),
              lower=lowbnd,upper=uppbnd)


best <- model$env$last.par.best
print(best)
rep <- sdreport(model)
print(summary(rep))

fit$objective
