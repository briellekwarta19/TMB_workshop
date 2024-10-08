#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
 {
  DATA_INTEGER(Nbatch);
  DATA_INTEGER(Nunit);
  DATA_IVECTOR(Batches);
  DATA_IVECTOR(Units);
  DATA_IVECTOR(Treat);
  DATA_VECTOR(Original);
  DATA_VECTOR(Final);
  DATA_INTEGER(IsRandom);                                                  // 0 is fixed only 1 is full model

  PARAMETER(Control);
  PARAMETER(Treatment);
  PARAMETER_VECTOR(EpsB);
  PARAMETER_VECTOR(EpsU);
  PARAMETER(logSigmaB);
  Type SigmaB = exp(logSigmaB);
  PARAMETER(logSigmaU);
  Type SigmaU = exp(logSigmaU);

  Type f;
  vector<Type>Prob(Nunit);

  f = 0;
  int u;
  int b;
  
  if(IsRandom == 1){ //fixed effect model
    for(u = 0; u< Nunit; u++){
      if(Treat(u) == 1){
        Prob(u) = 1.0/(1.0+exp(Control));
      }
      if(Treat(u) == 2){
       Prob(u) = 1.0/(1.0+exp(Treatment));
      }
      
      
    f-= dbinom(Final(u), Original(u), Prob(u), true);  //data,size,prob,true
    
    }

  }
  
  if(IsRandom == 2){ //random effect model
    for(u = 0; u< Nunit; u++){
      if(Treat(u) == 1){
        Prob(u) = 1.0/(1.0+exp(Control + EpsU(u) + EpsB(Batches(u)-1)));
      }
        
      if(Treat(u) == 2){
        Prob(u) = 1.0/(1.0+exp(Treatment + EpsU(u) + EpsB(Batches(u)-1)));
      }
      
      f-= dbinom(Final(u), Original(u), Prob(u), true);  //data,size,prob,true
      f-= dnorm(EpsU(u), Type(0), SigmaU, true);
    }
     
    for(b = 1; b < Nbatch; b++){
      f-= dnorm(EpsB(b), Type(0), SigmaB, true); 
    }
    
    }
  
  
  REPORT(f);
  ADREPORT(SigmaB);
  ADREPORT(SigmaU);
  REPORT(Prob);

  return f;
}




