//error function for penalizing the distnace from a reference 3D point coordinates
// e = X - X_r;


#include "ceres/rotation.h"
#include <iostream>
#include <fstream>
#include <iomanip>



class Reference3DPointError{
	public:
	Reference3DPointError(double * X_r, double weight){
		this->X_r= new double[3];
		this->X_r[0] = X_r[0];
		this->X_r[1] = X_r[1];
		this->X_r[2] = X_r[2];
		w=weight;
	}

	template <typename T>
	bool operator()(const T* const X, T* residuals)const{
		residuals[0] = T(w)*(T(X_r[0]) - X[0]);
		residuals[1] = T(w)*(T(X_r[1]) - X[1]);
		residuals[2] = T(w)*(T(X_r[2]) - X[2]);
		return true;
	}

	private:
	//reference point
	double * X_r;
	//weight
	double w;


};

