//error function for penalizing the distnace from a reference 3D point coordinates
// e = X - X_r;


#include "ceres/rotation.h"
#include <iostream>
#include <fstream>
#include <iomanip>



class ReferenceCamera7pError{
	public:
	ReferenceCamera7pError(double * X_r, double weight){
		this->X_r= new double[3];
		this->X_r[0] = X_r[0];
		this->X_r[1] = X_r[1];
		this->X_r[2] = X_r[2];
		w=weight;
	}

	template <typename T>
	bool operator()(const T* const X, T* residuals)const{
		residuals[0] = T(w)*(T(X_r[0]) - X[4]);
		residuals[1] = T(w)*(T(X_r[1]) - X[5]);
		residuals[2] = T(w)*(T(X_r[2]) - X[6]);
		return true;
	}

	private:
	//reference point
	double * X_r;
	//weight
	double w;


};

class ReferenceCamera6pError{
	public:
	ReferenceCamera6pError(double * X_r, double weight){
		this->X_r= new double[3];
		this->X_r[0] = X_r[0];
		this->X_r[1] = X_r[1];
		this->X_r[2] = X_r[2];
		w=weight;
	}

	template <typename T>
	bool operator()(const T* const X, T* residuals)const{
		residuals[0] = T(w)*(T(X_r[0]) - X[3]);
		residuals[1] = T(w)*(T(X_r[1]) - X[4]);
		residuals[2] = T(w)*(T(X_r[2]) - X[5]);
		return true;
	}

	private:
	//reference point
	double * X_r;
	//weight
	double w;


};

