library(stats4)

# the negative log-likelihood 
minusLL <- function(logLinf,logK,t0,logSigma,Ages,Lengths)
 {
  # Extract the parameters
  Linf <- exp(logLinf); Kappa <- exp(logK); Sigma <- exp(logSigma)
  
  # make the model predictions
  Pred <- Linf*(1.0-exp(-Kappa*(Ages-t0)))
  
  # Compute the negative log-likelihood
  NegLogL <- -1*sum(dnorm(Lengths,Pred,Sigma,TRUE))
  
  return(NegLogL)

 }  

 par(mfrow=c(2,2))


 # Psuedo data (always good for testing)
 set.seed(781010)
 Nsample <- 20
 Ages <- sample(c(1:20),Nsample,replace=TRUE)
 Lengths <- 100*(1.0-exp(-0.2*(Ages-0.1)))+rnorm(Nsample,0,5)
 plot(Ages,Lengths,xlab="Ages",ylab="Lengths")
 
 # fit the model
 start <- list(logLinf=log(80),logK=log(0.15),t0=0,logSigma=1)
 fixed <- list(Ages=Ages,Lengths=Lengths)
 mleOutput <- mle(minusLL,start=start,fixed=fixed)
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
 BestLL <- -as.numeric(logLik(mleOutput))
 print(BestLL)

 # Now do a likelihood profile
 sigs <- seq(from=3,to=10,by=0.1)
 Nprof <- length(sigs)
 NegLogLiks <- rep(0,Nprof)
 for (Isigma in 1:Nprof){
   # fit the model but FIX logSigma
   start <- list(logLinf=log(80),logK=log(0.15),t0=0)
   fixed <- list(logSigma=log(sigs[Isigma]),Ages=Ages,Lengths=Lengths)
   mleOutput <- mle(minusLL,start=start,fixed=fixed)
   NegLogLiks[Isigma]<- -as.numeric(logLik(mleOutput))-BestLL
  } 
 print(NegLogLiks)
   
 plot(sigs,NegLogLiks,xlab="Sigma",ylab="Negative log-likelihood",type="l",lty=1,ylim=c(0,10))
 # Indicate 95% CI
 abline(h=3.84/2)#3.84 = chi squared with 1 degree of freedom 
 points(Sigma,0,pch=16)
 
 #if 2dim use chi quared 2 degrees of freedom
 
