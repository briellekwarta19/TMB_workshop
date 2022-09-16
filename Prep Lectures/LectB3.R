library(here)

Nsize <- 3; Narea <- 2;Ndim <- Nsize*Narea

# Read in the marices 
DataFile <- here::here("Prep Lectures", "LectB3.txt")
S <- matrix(scan(DataFile,skip=1,n=Ndim*Ndim,quiet=T),nrow=Ndim,ncol=Ndim,byrow=T)
A <- matrix(scan(DataFile,skip=Ndim+2,n=Ndim*Ndim,quiet=T),nrow=Ndim,ncol=Ndim,byrow=T)
X <- matrix(scan(DataFile,skip=2*Ndim+3,n=Ndim*Ndim,quiet=T),nrow=Ndim,ncol=Ndim,byrow=T)

# Identity matrix - will be used often
I <- matrix(0,ncol=Ndim,nrow=Ndim); diag(I) = 1

FindEqn <- function(FF=0)
 {
  
  # Specify the H matrix
  H <- matrix(0,ncol=Ndim,nrow=Ndim); diag(H) = 1
  H[Nsize,Nsize] <- 1.0 - 0.2/0.8*FF
  H[Nsize+Nsize,Nsize+Nsize] <- 1.0 - FF
  
  # Multiply the matrices
  Mat1 <- S
  Mat2 <- (X %*% (A %*% (S %*% (H %*% S))))

  # Matrix inversion multiplied by a recruitment vector with a 1 in the first row
  Neqn <- solve(I-Mat2)[,1] #solves for matrix inversion
  #print(Neqn)
  
  # Code to chck the equilibrium is correct
  #Test <- S %*% Neqn
  #Test <- H %*% Test
  #Test <- S %*% Test
  #Test <- A %*% Test
  #Test <- X %*% Test
  #Test[1] <- Test[1] + 1
  #print(Test-Neqn)
  
  # Compute the yield
  Catch <- sum((I-H) %*% (Mat1 %*% Neqn))
  
  # Return
  Outs <- NULL
  Outs$FF <- FF #fishing effort
  Outs$Neqn <- Neqn #equilibrium
  Outs$Catch <- Catch
  return(Outs)
 }  


# =================================================================================

# Part 1 (yield vs F)
FFs <- seq(from=0,to=1,by=0.1)
Yields <- rep(0,length(FFs))
Spawn <- rep(0,length(FFs))
for (II in 1:length(FFs))
 {
  ModelPrj <- FindEqn(FFs[II])
  Yields[II] <- ModelPrj$Catch
  Spawn[II] <- sum(ModelPrj$Neqn[c(3,6)])
 }
par(mfrow=c(2,1),oma=c(1,1,1,1))
plot(FFs,Yields,xlab="Fishing effort",ylab="Yield",type="l",lty=1)
plot(FFs,Spawn,xlab="Fishing effort",ylab="Spawners",type="l",lty=1)




