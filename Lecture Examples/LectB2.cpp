#include <TMB.hpp>


template<class Type>
Type posfun(Type x, Type eps, Type &pen)
{
  if ( x >= eps ){
    return x;
  } else {
    pen += Type(0.01) * pow(x-eps,2);
    return eps/(Type(2.0)-x/eps);
  }
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(C); //Vector of real numbers
  DATA_VECTOR(I); //data basic data
  int n = C.size(); //n = length of C

  PARAMETER(logR); //log spaced parameters 
  PARAMETER(logK);
  PARAMETER(logQ);
  PARAMETER(logSigma);
  PARAMETER_VECTOR(FF);
  Type r = exp(logR); // Type r informs that r is a real number
  Type k = exp(logK);
  Type q = exp(logQ);
  Type sigma = exp(logSigma);
  
  //Temporary variables 
  int n1 = 0; 
  n1 = n + 1;
  vector<Type> B(n1);
  vector<Type> Ihat(n);
  vector<Type> Chat(n);
  vector<Type> ExpOut(n);
  Type f; // variable f is the objective function - the negatie log likelihood
  B(0) = k; //k = carry capacity, setting first element of the array to k
  
  //equations:
  for(int t=0; t<n; t++)
  {
    //exploitation rate: logit function of ff, transforms real # range to 0-1
    Type Expl = 1.0/(1.0+exp(-FF(t))); 
    B(t+1) = B(t) + r*B(t)*(1-B(t)/k) - Expl*B(t); //Biomass
    Chat(t) = Expl*B(t); //predicted catch
    ExpOut(t) = Expl; //exploitation rate
    Ihat(t) = q*B(t); //predicted CPUE
  }
  f = -sum(dnorm(log(C), log(Chat), Type(0.05), true)); //negative log likelihood
  f -= sum(dnorm(log(I), log(Ihat), sigma, true)); 
  
  //above is equivalent to f = f-sum(.....)
   
  //f = -sum(dlnorm(C, log(Chat), Type(0.05), true));
  //f -= sum(dlnorm(I, log(Ihat), sigma, true));

  ADREPORT(log(B)); // uncertainty -printing uncertainties 
  REPORT(f);        // plot
  REPORT(B);        // plot
  REPORT(Chat);        // plot
  REPORT(ExpOut);        // plot
  REPORT(Ihat);     // plot

  return f;
}




// dvariable posfun(const dvariable&x,const double eps,dvariable& pen)
//     {
//       if (x>=eps) {
//         return x;
//       } else {
//         pen+=.01*square(x-eps);
//         return eps/(2-x/eps);
// } }


