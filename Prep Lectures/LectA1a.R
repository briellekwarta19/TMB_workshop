library(stats4)

Ages=NULL
Lengths=NULL

# the negative log-likelihood for the growth function in lecture 1
minusLL <- function(logLinf,logK,t0,logSigma)
 {
  # Extract the parameters, log space to keep values positive
  Linf <- exp(logLinf); Kappa <- exp(logK); Sigma <- exp(logSigma)
  
  # make the model predictions
  Pred <- Linf*(1.0-exp(-Kappa*(Ages-t0)))
  
  # Compute the negative log-likelihood
  NegLogL <- -1*sum(dnorm(Lengths,Pred,Sigma,TRUE))
  
  return(NegLogL)

 }  

 # Psuedo data (always good for testing)
 set.seed(781010)
 Nsample <- 5000
 #generating data
 Ages <<- sample(c(1:20),Nsample,replace=TRUE) 
 Lengths <<- 100*(1.0-exp(-0.2*(Ages-0.1)))+rnorm(Nsample,0,5)
 plot(Ages,Lengths,xlab="Ages",ylab="Lengths")
 
 # fit the model
 #start = a list of parameters to be estimated and initial values
 start <- list(logLinf=log(80),logK=log(0.15),t0=0,logSigma=1) 
 fixed <- NULL #parameters we wish to fix
 fixed <- list(Ages = Ages, Lengths = Lengths) #parameters we wish to fix
 mleOutput <- mle(minusLL,start=start) 
 print(summary(mleOutput))
 
 # Extract the parameters
 Linf <- exp(coef(mleOutput)[1])
 Kappa <- exp(coef(mleOutput)[2])
 t0 <- coef(mleOutput)[3]
 Sigma <- exp(coef(mleOutput)[4])
 cat("parameter estimates\n")
 cat("Linf = ",Linf,"\n")
 cat("Kappa = ",Kappa,"\n")
 cat("t0 = ",t0,"\n")
 cat("Sigma = ",Sigma,"\n")
 
