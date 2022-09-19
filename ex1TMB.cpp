#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(Age);
  DATA_VECTOR(Length);
  DATA_INTEGER(Model);
  
  int n = Length.size();

  PARAMETER(logLinf);
  PARAMETER(loga50);
  PARAMETER(a0);
  PARAMETER(logk);
  PARAMETER(logSigma);
  
  vector<Type> Lafit(n);

  Type neglogL = 0.0;
  Type Linf = exp(logLinf); 
  Type a50 = exp(loga50); 
  Type k = exp(logk);

  Type delta = (Linf*0.95) - a50;
  REPORT(delta);
  
  if(Model = 1){
    Lafit = Linf * 1/(1 + exp(-log(19)*(Age- a50)/delta));
  }

  else{
    Lafit = Linf * (1-exp(-k*(Age-a0)));
  }

  neglogL = -sum(dnorm(Length, Lafit, exp(logSigma), true));
  
  return neglogL;
  
}
