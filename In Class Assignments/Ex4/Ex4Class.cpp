#include <TMB.hpp>

template <class Type> Type square(Type x){return x*x;}

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_INTEGER(Nyear)
  DATA_INTEGER(Nclass)
  DATA_VECTOR(Length)
  DATA_VECTOR(Weight)
  DATA_MATRIX(X);
  DATA_VECTOR(S)
    DATA_VECTOR(SurveyS)
    DATA_SCALAR(M)
    DATA_VECTOR(CWObs)
    DATA_MATRIX(CALObs);
  DATA_SCALAR(Neff)
    DATA_VECTOR(BioIndex)
    DATA_SCALAR(BioSig)
    DATA_INTEGER(Nproj)
    DATA_SCALAR(Fproj)
    
    // End of data section
    
    PARAMETER(dummy);
  PARAMETER(LogRbar);
  PARAMETER_VECTOR(LogNinit);
  PARAMETER_VECTOR(LogFullF);
  PARAMETER_VECTOR(Eps);
  
  matrix<Type> N(Nyear+Nproj+1,Nclass);
  matrix<Type> F(Nyear+Nproj,Nclass);
  matrix<Type> Z(Nyear+Nproj,Nclass);
  matrix<Type> CAL(Nyear+Nproj,Nclass);
  vector<Type> CW(Nyear+Nproj);
  vector<Type> BioPred(Nyear+Nproj);
  
  Type CALtot;
  
  Type Penal;
  Type LikeCatch;
  Type LikeBio;
  Type LikeCAL;
  Type obj_fun;
  
  // End of specifications section
  // =============================
  
  // First set F and Z by size-classs (note that Fproj applies after year Nyear)
  for (int Iyear=0; Iyear<Nyear+Nproj; Iyear++)
    for (int Iclass=0;Iclass<Nclass;Iclass++)
    {
      if (Iyear < Nyear)
        F(Iyear,Iclass) = exp(LogFullF(Iyear))*S(Iclass);
      else
        F(Iyear,Iclass) = Fproj*S(Iclass);
      Z(Iyear,Iclass) = M + F(Iyear,Iclass);
    }
    
    // Now set the N matrix
    for (int Iclass=0;Iclass<Nclass;Iclass++) N(0,Iclass) = exp(LogNinit(Iclass));
  for (int Iyear=0;Iyear<Nyear+Nproj;Iyear++)
  {
    // Catch-at-length
    CALtot = 0; CW(Iyear) = 0;
    for (int Iclass=0;Iclass<Nclass;Iclass++)
    {
      CAL(Iyear,Iclass) = F(Iyear,Iclass)/Z(Iyear,Iclass)*N(Iyear,Iclass)*(1.0-exp(-Z(Iyear,Iclass)));
      CALtot += CAL(Iyear,Iclass);
      CW(Iyear) += Weight(Iclass)*CAL(Iyear,Iclass);
    }
    CALtot = (CAL.row(Iyear)).sum();
    for (int Iclass=0;Iclass<Nclass;Iclass++) CAL(Iyear,Iclass) /= CALtot;
    
    // Numbers-at-age
    for (int Iclass=0;Iclass<Nclass;Iclass++)
    {
      N(Iyear+1,Iclass) = 0;
      for (int Jclass=0;Jclass<Nclass;Jclass++)
        N(Iyear+1,Iclass) += N(Iyear,Jclass)*exp(-Z(Iyear,Jclass))*X(Iclass,Jclass);
    }
    
    // Recruitment (watch for the index for Eps - and N)
    N(Iyear+1,0) += exp(LogRbar)*exp(Eps[Iyear]);
  }
  
  
  //Catch Likelihood, 0.05 = s.d. of catch weight
  Type SS = 0;
  // for (int Iyear=0; Iyear<Nyear+Nproj; Iyear++)
  // {
  //   SS += pow((log(CW(Iyear)) - log(CWObs(Iyear))), Type(2.0));
  // }
  // 
  // LikeCatch = (1.0 / (2.0 * square(0.05))) * SS;
  
  // Biomass predictions
  for (int t=0;t<Nyear;t++){
    BioPred(t) = 0.0;
    for (int i=0;i<Nclass;i++){
      BioPred(t) += SurveyS(i)*Weight(i)*N(t,i);
    }
  }
  
  // Index Likelihood (BioIndex)
  SS = 0.0;
  Type q = 0.0;
  
  for (int t=0;t<Nyear;t++){
    q +=log(BioIndex(t)/BioPred(t));
  }
  
  q = exp(q/25.0);
  
  for (int t=0;t<Nyear;t++){
    SS += (log(BioIndex(t)) - log(q*BioPred(t)));
  }
  
  LikeBio = 1/(2*pow(0.05,2))*SS;
  
  // CAL Likelihood
  LikeCAL = 0;
  
  // Recruitment penalty (include years after Nyear)
  Penal = 0;
  
  obj_fun = dummy*dummy + LikeCatch+LikeBio+LikeCAL+Penal;
  
  // Stuff to report
  REPORT(N);
  REPORT(Z);
  REPORT(LikeCatch);
  REPORT(LikeBio);
  REPORT(LikeCAL);
  REPORT(Penal);
  REPORT(BioPred);
  REPORT(obj_fun);
  
  return(obj_fun);
}
