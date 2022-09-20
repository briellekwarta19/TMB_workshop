#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
	DATA_VECTOR(C); //citations
  DATA_INTEGER(Model); // Model type

	int t = C.size(); //number of years

	//define parameters
	PARAMETER(loga); // scalar
	PARAMETER(logb); // scalar
	PARAMETER(logc); // scalar
	PARAMETER(logd); // scalar
	PARAMETER(loge); // scalar

	Type a = exp(loga);
	Type b = exp(logb);
	Type c = exp(logc);
	Type d = exp(logd);
	Type e = exp(loge);

	//Predictions vector
	vector<Type> Chat(t);
	
	//define objective function
	Type nll = 0.0;

	//Calculate citation predictions
	for (int i = 0; i < t; i++){
	  if(Model == 1){
	    Chat(i) = a*i;
	  }
    if(Model == 2){ // CRAZY MODEL 
      Chat(i)=a*sin(b)+ c* pow(i,2) + d*i + e;
    }
    
    if(Model == 3){ // m shape curve
      Chat(i)=a*(i+b)*pow((i+c),2)*pow((i+d),3) + e;
    }
    
    if(Model == 4){ // power
      Chat(i)=a*(i+b)*pow((i+c),2) + d;
    }
    
    if(Model == 5){ // sin
     Chat(i)=a*(i+b)+ sin(i*c);
    }
    
    
	 }
	 
  nll = -sum(dpois(C, Chat, true));

	return(nll);
}
