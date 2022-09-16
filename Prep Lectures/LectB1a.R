# NOTE: Something is up with code, couldn't get mle to work with minusLL function initially
# Compared with code from past workshop, seems to be related to including C and I are 
# in the function/fixed list (weren't in the code from last year)
# Also, the Type=1 wasn't included in the set of fixed parametrs
# Below is the function definition for minusLL in from the last workshop
# minusLL <- function(logK,logr, Type=1)

library(stats4)

# the negative log-likelihood 
minusLL <- function(logK, logr, Type = 1)
{
  # Extract the parameters
  K <- exp(logK); r <- exp(logr); Nyear <- length(C); B <- rep(0, Nyear+1)
  
  # Calculate the biomass series
  B[1] <- K # initial biomass = K, meaning unfished
  
  for (Iyear in 1:Nyear){
    B[Iyear+1] <- B[Iyear] + r*B[Iyear]*(1.0-B[Iyear]/K) - C[Iyear]
    # set lower limit for biomass
    if (B[Iyear] < 0.01) B[Iyear] <- 0.01
  }
  
  print(B)
  # ML estimate for q
  qhat <- 0; bot <- 0;
  for (Iyear in 1:Nyear){
    
    qhat <- qhat + log(I[Iyear]/B[Iyear])
    bot <- bot + 1
    
    }
  qhat <- exp(qhat/bot)
  print(qhat)
  # Compute the sum of squared errors and hence the ML estimate for Sigma
  SS <-  0;
  for (Iyear in 1:Nyear){
    
    SS <- SS + (log(I[Iyear]) - log(qhat*B[Iyear]))^2
    
  }
  Sigma <- sqrt(SS/Nyear)
  print(Sigma)
  # Concentrated likelihood
  NegLogLike <- Nyear*log(Sigma);
 
  # material to return
  Outs <- NULL
  Outs$B <- B
  Outs$CpueHat <-B*qhat
  Outs$r <- r
  Outs$K <- K
  print(Type)
  if (Type == 1) return(NegLogLike)
  if (Type == 2) return(Outs)
 }  

# the negative log-likelihood 
minusLL2 <- function(logMSY,logr,  Type=1)
{
  # Extract the parameters
  MSY <- exp(logMSY); r <- exp(logr); K <- 4*MSY/r; Nyear <- length(C); B <- rep(0, Nyear+1)

  # biomass series
  B[1] <- K
  for (Iyear in 1:Nyear)
  {
    B[Iyear+1] <- B[Iyear] + r*B[Iyear]*(1.0-B[Iyear]/K) - C[Iyear]
    if (B[Iyear] < 0.01) B[Iyear] <- 0.01
  }  
  
  # ML estimate for q
  qhat <- 0; bot <- 0;
  for (Iyear in 1:Nyear){ 
    
    qhat <- qhat + log(I[Iyear]/B[Iyear]); bot <- bot + 1 
    
    }  
  qhat <- exp(qhat/bot)
  
  # Concentrated likelihood
  SS <-  0;
  for (Iyear in 1:Nyear){
    
    SS <- SS + (log(I[Iyear]) - log(qhat*B[Iyear]))^2
    
  }
  Sigma <- sqrt(SS/Nyear)
  
  NegLogLike <- Nyear*log(Sigma);
  #cat(r,K,NegLogLike,"\n")
  
  # material to return
  Outs <- NULL
  Outs$B <- B
  Outs$CpueHat <-B*qhat
  Outs$r <- r
  Outs$K <- K
  
  if (Type==1) return(NegLogLike)
  if (Type==2) return(Outs)
}  

par(mfrow=c(2,2))


# Catch and effort data
C <- c(15.9,25.7,28.5,23.7,25.0,33.3,28.2,19.7,17.5,19.3,21.6,23.1,22.5,22.5,23.6,29.1,14.4,13.2,28.4,34.6,37.5,25.9,25.3) 
I <- c(61.89,78.98,55.59,44.61,56.89,38.27,33.84,36.13,41.95,36.63,36.33,38.82,34.32,37.64,34.01,32.16,26.88,36.61,30.07,30.75,23.36,22.36,21.91)
Yr <- seq(from=1967,to=1989)


# Part 1 (fit the model and plot diagnostics)
# ===========================================

# fit the model
start <- list(logK=log(400), logr=log(0.3))
fixed <- list( Type = 1)
mleOutput <- mle(minusLL,start=start,fixed=fixed)
print(summary(mleOutput))

# Extract the model parameters
logK <- coef(mleOutput)[1]
logr <- coef(mleOutput)[2]
Output <- minusLL(logK,logr,Type=2)

# Plot the fit
plot(c(Yr,1990),Output$B,ylim=c(0,max(Output$B)*1.05),type="l",xlab="Year",ylab="Biomass ('000t)")
plot(c(Yr,1990),Output$CpueHat,ylim=c(0,max(I)*1.05),type="l",xlab="Year",ylab="Biomass ('000t)")
points(Yr,I,pch=16)

# Part 2 (profile likelihood for MSY)
# ===================================

# fit the model
start <- list(logMSY=log(400*0.3/2),logr=log(0.01))
fixed <- list(Type = 1)
mleOutput <- mle(minusLL2,start=start,fixed=fixed)
print(summary(mleOutput))


#start <- list(logMSY=log(400*0.3/2),logr=log(0.01))
#fixed <- list(C=C,I=I)
#mleOutput <- mle(minusLL2,start=start,fixed=fixed)
#print(summary(mleOutput))

