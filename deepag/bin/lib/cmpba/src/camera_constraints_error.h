// Camera constraints error
// for any number of parameters
// number of parameters specified by ncp
// block - [fx,x0,y0,ar,r0,r1,q1,q2,q3,t0,t1,t2]
class CameraConstraintsError{
	public:
	CameraConstraintsError(double * constr, double * weights, int ncp){
		this->ncp = ncp;
		this->constr = constr;
		this->weights = weights;

	}
	

	template <typename T>
	bool operator()(const T* const p,  T* residuals)const{
		for(int i=0;i<ncp;i++){
			residuals[i]=(constr[i]-p[i])*weights[i];		
		}
		return true;
	}

	private:
	double * constr;
	double * weights;
	int ncp;
};
