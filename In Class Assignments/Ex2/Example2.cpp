#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(NC); // Predator vector
  DATA_MATRIX(NP); //Prey counts matrix 
  DATA_MATRIX(PY); //Consumption counts matrix
  DATA_INTEGER(Model); // Model 
  
  int n = NC.size(); //n = length of C
  
  PARAMETER_VECTOR(loga);
  //PARAMETER_VECTOR(logB)
  //PARAMETER_VECTOR(loggamma);

  vector<Type> a = exp(loga); // Type r informs that r is a real number
  //vector <Type> B = exp(logB);
  //vector <Type> gamma = exp(loggamma);
  
  //Temporary variables 
  
  matrix<Type> Phat(n,3);
  Type ss = 0.0; //sum of squares

  //equations:
  //MODEL 1
  if(Model == 1){
    for(int i= 0; i < n; ++i){
      for(int j= 0; j < 2; ++j){
        Phat(i,j) = a(i)*NC(j);
      }
    }
  }
  
  //MODEL 2
  // if(Model == 2){
  //   for(int i= 0; i < n; ++i){
  //     for(int j= 0; j < 2; ++j){
  //       Phat(i,j) = (a(i)*NC(j))/ (1+ B(i)*NP(i,j));
  //     }
  //   }
  // }
  
  //MODEL 3
  // if(Model == 3){
  //   for(int i= 0; i < n; ++i){
  //     for(int j= 0; j < 2; ++j){
  //       Phat(i,j) = (a(i)*NC(j)*pow(NP(i,j), gamma - 1)) / (1+ B(i)*pow(NP(i,j), gamma));
  //     }
  //   }
  // }
  // 
  
  for(int i = 0; i < 2; ++i){
    for(int j = 0; j < 2 ; ++j){
      
      ss += pow(log(PY(i,j)) - log(Phat(i,j)), Type(2)); 
    }
  }
  

  return -ss;
}
