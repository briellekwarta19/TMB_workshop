#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
	DATA_VECTOR(NC);
	DATA_MATRIX(NP);
	DATA_MATRIX(PY);
	DATA_INTEGER(Model);
	int n = NC.size();

	//define parameters
	PARAMETER_VECTOR(loga);

	vector<Type> a = exp(loga);

	//define matrix for predictions
	
	matrix<Type> Phat(n,3);
	
	//define objective function
	Type ss = 0.0; //sum of squares

	//loop to calculate preds
	for (int i = 0; i < n; ++i){
	  for(int j = 0; j < 2; ++j){
	 		Phat(i,j) = a(j)*NC(i);
	 	}
	 }

	//double loop to create sum of squares
	for(int i = 0; i < n; ++i){
		for(int j = 0; j < 2; ++j){
			ss += pow(log(PY(i,j))-log(Phat(i,j)),Type(2));
		}
	}

	return(ss);
}
