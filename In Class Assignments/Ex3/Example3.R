setwd("In Class Assignments/Ex3")

library(tidyverse)
library(gridExtra)

data <- read.table("ex3.dat", header=FALSE)
colnames(data) <- c("Year", "Google", "Web")

p1 <- ggplot(data) +
  geom_path(mapping = aes(x = Year, y = Google))

p2 <- ggplot(data)+
  geom_path(mapping = aes(x = Year, y = Web))
  
grid.arrange(p1,p2)

require(TMB)

compile("Example3.cpp", flags="-Wno-ignored-attributes")
dyn.load(dynlib("Example3"))

################################################################################
#### Model 1 ####
data_mod1 <- list(C = as.vector(data$Google), Model = 1)

parameter1 <- list(loga= 1, logb = 1, logc = 1, logd = 1, loge = 1)

map1 <- list(logb=factor(NA),logc=factor(NA), logd=factor(NA), loge = factor(NA))

model <- MakeADFun(data_mod1, 
                   map = map1,
                   parameter1,
                   DLL="Example3")

fit <- nlminb(model$par, model$fn, model$gr)
nll <- print(fit$objective)
best <- model$env$last.par.best
rep <- sdreport(model)
print(summary(rep))

k <- 1
n <- length(data$Year)
AIC <- 2*k*n/ (n-k-1) + 2*nll

#### Model 2 ####
data_mod2 <- list(C = as.vector(data$Google), Model = 2)

parameter2 <- list(loga= 1, logb = 1, logc = 1, logd = 1, loge = 1)

model <- MakeADFun(data_mod2, 
                   parameter2,
                   DLL="Example3")

fit <- nlminb(model$par, model$fn, model$gr)
nll <- (fit$objective)
best <- model$env$last.par.best
rep <- sdreport(model)
print(summary(rep))


k <- 5
n <- length(data$Year)
AIC <- 2*k*n/ (n-k-1) + 2*nll
AIC

#### Model 3 ####
data_mod3 <- list(C = as.vector(data$Google), Model = 3)

parameter3 <- list(loga= 1, logb = 1, logc = 1, logd = 1, loge = 1)

model <- MakeADFun(data_mod3, 
                   parameter3,
                   DLL="Example3")

fit <- nlminb(model$par, model$fn, model$gr)
nll <- (fit$objective)
best <- model$env$last.par.best
rep <- sdreport(model)
print(summary(rep))


k <- 5
n <- length(data$Year)
AIC <- 2*k*n/ (n-k-1) + 2*nll
AIC

#### Model 4 ####
data_mod4 <- list(C = as.vector(data$Google), Model = 4)

parameter4 <- list(loga= 1, logb = 1, logc = 1, logd = 1, loge = 1)
map4 <- list(loge=factor(NA))

model <- MakeADFun(data_mod4, 
                   parameter4,
                   map= map4,
                   DLL="Example3")

fit <- nlminb(model$par, model$fn, model$gr)
nll <- (fit$objective)
best <- model$env$last.par.best
rep <- sdreport(model)
print(summary(rep))


k <- 4
n <- length(data$Year)
AIC <- 2*k*n/ (n-k-1) + 2*nll
AIC


#### Model 5 ####
data_mod5 <- list(C = as.vector(data$Google), Model = 5)

parameter5 <- list(loga= 1, logb = 1, logc = 1, logd = 1, loge = 1)
map5 <- list(logd = factor(NA), loge=factor(NA))

model <- MakeADFun(data_mod5, 
                   parameter5,
                   map= map5,
                   DLL="Example3")

fit <- nlminb(model$par, model$fn, model$gr)
nll <- (fit$objective)
best <- model$env$last.par.best
rep <- sdreport(model)
print(summary(rep))

k <- 3
n <- length(data$Year)
AIC <- 2*k*n/ (n-k-1) + 2*nll
AIC

mod5_est <- rep(NA, n)
sum5 <- summary(rep)
a <- exp(sum5[1,1])
b <- exp(sum5[2,1])
c <- exp(sum5[3,1])

for(t in 1:n){
  mod5_est[t]=a*(t+b)+ sin(t*c)
}

data5 <- data

data5 <- cbind(data,mod5_est)

ggplot(data5) +
  geom_path(mapping = aes(x = Year, y = Google))+
  geom_path(mapping = aes(x = Year, y = mod5_est), col = 'red')


#### Model 5 2 data####
data_mod5 <- list(C = as.vector(data$Web), Model = 5)

parameter5 <- list(loga= 1, logb = 1, logc = 1, logd = 1, loge = 1)
map5 <- list(logd = factor(NA), loge=factor(NA))

model <- MakeADFun(data_mod5, 
                   parameter5,
                   map= map5,
                   DLL="Example3")

fit <- nlminb(model$par, model$fn, model$gr)
nll <- (fit$objective)
best <- model$env$last.par.best
rep <- sdreport(model)
print(summary(rep))

k <- 3
n <- length(data$Year)
AIC <- 2*k*n/ (n-k-1) + 2*nll
AIC

mod5_est <- rep(NA, n)
sum5 <- summary(rep)
a <- exp(sum5[1,1])
b <- exp(sum5[2,1])
c <- exp(sum5[3,1])

for(t in 1:n){
  mod5_est[t]=a*(t+b)+ sin(t*c)
}

data5 <- data

data5 <- cbind(data,mod5_est)

ggplot(data5) +
  geom_path(mapping = aes(x = Year, y = Web, col = 'WebofSci'))+
  geom_path(mapping = aes(x = Year, y = mod5_est, col = 'Model5'))




