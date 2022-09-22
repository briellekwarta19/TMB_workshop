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
  Type g;
  Type DevRes;
  Type ParEffect;
  vector<Type>Prob(Nunit);

  f = 0; g = 0;
  for (int Index=0; Index<Nunit; Index++)
   {
	  if (Treat(Index)==1) ParEffect = Control;
	  if (Treat(Index)==2) ParEffect = Control+Treatment;
	  ParEffect = ParEffect + EpsB(Batches(Index)-1) + EpsU(Index);
	  Prob(Index) = 1.0/(1.0+exp(ParEffect));
	
	  // Negative log-likelihood (estimated p)
	  DevRes = dbinom(Final(Index),Original(Index),Prob(Index),true);
	  // Negative log-likelihood (saturated p)
	  DevRes -= dbinom(Final(Index),Original(Index),Final(Index)/Original(Index),true);
	  f -= DevRes;
	  g -= DevRes;
  }
  // Deviance is 2 * (log(D|p_sat) - log(D|p_sat))
  Type Deviance = 2*f;
  // if (IsRandom==1) f -= sum(dnorm(EpsB,Type(0.0),SigmaB,true));
  if (IsRandom==1) 
    for (int j=0;j<Nbatch;j++) f -= dnorm(EpsB(j),Type(0.0),SigmaB,true);
  // if (IsRandom==1) f -= sum(dnorm(EpsU,Type(0.0),SigmaU,true));
  if (IsRandom==1) 
    for (int j=0;j<Nunit;j++) f -= dnorm(EpsU(j),Type(0.0),SigmaU,true);

  REPORT(g);
  REPORT(f);
  REPORT(Deviance);
  ADREPORT(SigmaB);
  ADREPORT(SigmaU);
  REPORT(Prob);

  return f;
}




