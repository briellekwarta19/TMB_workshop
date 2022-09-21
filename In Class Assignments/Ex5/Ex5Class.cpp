#include <TMB.hpp>

template <class Type> Type square(Type x){return x*x;}

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_INTEGER(m)
  DATA_IVECTOR(TT)
  DATA_INTEGER(Tmax)
  DATA_MATRIX(B)
  DATA_MATRIX(R)
  DATA_VECTOR(Phi0)

  // End of data section

  PARAMETER(dummy);
  PARAMETER(mu);
  PARAMETER(log_tau);
  PARAMETER_VECTOR(B0);
  PARAMETER_VECTOR(log_sigR);
  PARAMETER_VECTOR(eta);
  // End of the estimated parameters

  vector<Type> R0(m);
  vector<Type> SigR(m);
  vector<Type> h(m);
  Type tau;
  Type beta;
  Type BH;
  Type obj_fun;
  matrix<Type> BHout(Tmax,m);
  // End of the temporary variables

  // Transform the parameters
  R0 = B0/Phi0;
  tau = exp(log_tau);
  SigR = exp(log_sigR);

  obj_fun = 0;
 for (int k=0;k<m;k++) {

   // Extract beta and define h
   beta = mu + tau*eta(k);
   
   h(k) = (exp(beta + 0.2))/ (exp(beta + 1.0));
  
   for(int t = 0; t < Tmax; t++){
     BH = (4.0 * R0(k) * h(k) * B(k,t))/((Phi0(k)*R0(k)*(1.0-h(k)))+ (5.0*h(k)-1.0)*B(k,t));
      
     BHout(t,k) = BH;
   // Likelihood
   
   obj_fun += square(log(R(k,t))-log(BH) + square(SigR(k))/2.0)/(square(SigR(k))) + log(SigR(k));
   
   }
   
   obj_fun -= dnorm(eta(k), Type(0), Type(1), true);

   
  }



  obj_fun += dummy*dummy;
  ADREPORT(h);
  ADREPORT(R0);
  ADREPORT(tau);
  REPORT(BHout);

  return(obj_fun);
}
