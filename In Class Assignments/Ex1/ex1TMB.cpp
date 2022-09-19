#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(Age);
  DATA_VECTOR(Length);
  Data_Integer(Model);
  int n = Length.size();

  PARAMETER(Linf);
  PARAMETER(a50);
  PARAMETER(a0);
  PARAMETER(k);
  parameter(logsigma);
  
  vector<Type> Lafit(n);

  Type neglogL = 0.0;
  Type Linf = log(Linf); 
  Type a50 = log(a50); 
  Type k = log(k);

  Type delta = (Linf*0.95) - a50;
  
  if(Model = 1){
    Lafit = Linf * pow(1 + exp(-log(19)*(Age- a50)/delta),-1);
  }

  if(Model = 2){
    Lafit = Linf * (1-exp(-k*(Age-a0)));
  }

  neglogL = -sum(dnorm(Length, Lafit, exp(logSigma), true));
  
  return neglogL;
  
}
