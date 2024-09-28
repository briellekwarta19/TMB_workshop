library(TMB)

#Compilation. Returns 0 if compilation was successful
#setwd("bacri-sup/code")
TMB::compile("bacri-sup/code/linreg.cpp")

#load dynamic library (loading information from the model?)
dyn.load(dynlib("bacri-sup/code/linreg"))

#set seed
set.seed(123)

#say the linear regression model has a mean of a + bx
data <- list(y = rnorm(20) + 1:20, x = 1:20)

parameters <- list(a = 0, b = 0, tsigma = 0)

#instruct TMB to create likelihood function
obj_linreg <- MakeADFun(data, parameters, DLL = "linreg", silent = TRUE)

#Optimization of objective function with nlminb (negative log likelihood minimizer)
mod_linreg <- nlminb(start = obj_linreg$par, #starting initial parameters
                    objective = obj_linreg$fn, #objective function to minimize
                    gradient = obj_linreg$gr, #gradient: tells us which direction to step in likelihood estimation (derivative of objective in respect to each parameter)
                    hessian = obj_linreg$he) #hessian: tells us how far to step in the likelihood estimation (second derivative of objective in respect to each parameter)

mod_linreg$par

#returning maximum likelihood estimates and standard errors of parameters in terms of which nll is parameterized
sdreport(obj_linreg, par.fixed = mod_linreg$par)

#standard error reports
summary(sdreport(obj_linreg, par.fixed = mod_linreg$par), select = "report")

#displaying estimation results from lm function for comparison

rbind(
  "lm" = lm(y~x, data = data)$coeff, #linear regression in R
  "TMB" = mod_linreg$par[1:2] #intercept and slope from TMB fit
)

