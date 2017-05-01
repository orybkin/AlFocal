//an implementation of a reprojection error function according to Lourakis
//used parameters are k,q,t
//k = fx,x0

#include "ceres/rotation.h"
#ifdef MEXDEBUG
#include <mex.h>
#endif

using namespace ceres;

// Reprojection functions based on Lourakis'SBA
// |x|   |fx   s   x0|
// |y| = |0  ar*fx y0|*R*[1|t]*X
// |z|   |0    0    1|
//
// |u|   |x/z|
// |v| = |y/z|
//
// Camera has 12 pramaeters or 14 with radial distortion
// only 11 and 13 are adjusted at max, because only vector part of quaternion is adjusted, the scalar component is 
// calculated as q0 = sqrt(q1^2+q2^2+q3^2)

// 0 - parameters adjusted or fixed based on parameter mask
// 13 parameters
// blocks - [q1,q2,q3][t0,t1,t2][fx][ar][x0,y0][s][r1,r2]
// all blocks are fixable through parmask
class ReprojectionErrorCustom{
	public:
	ReprojectionErrorCustom(double observed_x, double observed_y, double * q0, double weight)
	:observed_x(observed_x),observed_y(observed_y){
		this->q0[0] = q0[0];
		this->q0[1] = q0[1];
		this->q0[2] = q0[2];
		this->q0[3] = q0[3];
		this->weight = weight;
	}
	

	template <typename T>
	bool operator()(const T* q, const T* C, const T* fx, const T* xy0, const T* ar, const T* s, const T* r, const T* const X, T* residuals)const{
		T w = ceres::sqrt(T(1)-q[0]*q[0]-q[1]*q[1]-q[2]*q[2]);
		const T qlocal[4]= {w,q[0],q[1],q[2]};
		T qtotal[4];
		const T tq0[4] = {T(q0[0]), T(q0[1]), T(q0[2]),T(q0[3])};
		QuaternionProduct( qlocal, tq0, qtotal);
		T translatedX[3];
		translatedX[0] = X[0] - C[0];
		translatedX[1] = X[1] - C[1];
		translatedX[2] = X[2] - C[2];
		T rotatedX[3];
		QuaternionRotatePoint(qtotal, translatedX, rotatedX);
		T x,y;
		x = rotatedX[0] / rotatedX[2];
		y = rotatedX[1] / rotatedX[2];
		T p_norm2 = x*x+y*y;
		T dist = T(1)+r[0]*p_norm2+r[1]*p_norm2*p_norm2;
		x=dist*x;
		y=dist*y;
		x = fx[0]*x+s[0]*y+xy0[0];
		y = ar[0]*fx[0]*y+xy0[1];
		residuals[0] = T(weight)*(T(observed_x) - x);
		residuals[1] = T(weight)*(T(observed_y) - y);
		return true;
	}

	private:
	double observed_x;
	double observed_y;
    double q0[4];
	double weight;
};

// 1 - all parameters adjusted
// 11 parameters
// block - [q1,q2,q3,t0,t1,t2,fx,x0,y0,ar,s]
class ReprojectionErrorKQT{
	public:
		ReprojectionErrorKQT(double observed_x, double observed_y, double * q0, double weight)
	:observed_x(observed_x),observed_y(observed_y){
		this->q0[0] = q0[0];
		this->q0[1] = q0[1];
		this->q0[2] = q0[2];
		this->q0[3] = q0[3];
		this->weight = weight;
	}
	

	template <typename T>
	bool operator()(const T* kqt, const T* const X, T* residuals)const{
		const T * k = kqt+6;
		const T * q = kqt;
		const T * t = kqt+3;
		T w = ceres::sqrt(T(1)-q[0]*q[0]-q[1]*q[1]-q[2]*q[2]);
		const T qlocal[4]= {w,q[0],q[1],q[2]};
		T qtotal[4];
		const T tq0[4] = {T(q0[0]), T(q0[1]), T(q0[2]),T(q0[3])};
		QuaternionProduct( qlocal, tq0, qtotal);
		T rotatedX[3];
		QuaternionRotatePoint(qtotal, X, rotatedX);
		T translatedX[3];
		translatedX[0] = rotatedX[0]+t[0];
		translatedX[1] = rotatedX[1]+t[1];
		translatedX[2] = rotatedX[2]+t[2];
		T x,y;
		x = k[0]*translatedX[0]+k[4]*translatedX[1]+k[1]*translatedX[2];
		y = k[3]*k[0]*translatedX[1]+k[2]*translatedX[2];
		residuals[0] = T(weight)*(T(observed_x) - x / translatedX[2]);
		residuals[1] = T(weight)*(T(observed_y) - y / translatedX[2]);
		return true;
	}

	private:
	double observed_x;
	double observed_y;
    double q0[4];
	double weight;
};

// 2 - K is fixed, provided calibration is used and kept constant
// 6 parameters
// block - [q1,q2,q3,t0,t1,t2]
class ReprojectionErrorQT{
	public:
		ReprojectionErrorQT(double observed_x, double observed_y, double * q0, double * k, double weight)
	:observed_x(observed_x),observed_y(observed_y){
		this->q0[0] = q0[0];
		this->q0[1] = q0[1];
		this->q0[2] = q0[2];
		this->q0[3] = q0[3];
		this->k[0] = k[0];
		this->k[1] = k[1];
		this->k[2] = k[2];
		this->k[3] = k[3];
		this->k[4] = k[4];
		this->weight = weight;
	}
	

	template <typename T>
	bool operator()(const T* const qt, const T* const X, T* residuals)const{
		const T * q = qt;
		const T * t = qt+3;
		T w = ceres::sqrt(T(1)-q[0]*q[0]-q[1]*q[1]-q[2]*q[2]);
		const T qlocal[4]= {w,q[0],q[1],q[2]};
		T qtotal[4];
		const T tq0[4] = {T(q0[0]), T(q0[1]), T(q0[2]),T(q0[3])};
		QuaternionProduct(qlocal, tq0, qtotal);
		T rotatedX[3];
		QuaternionRotatePoint(qtotal, X, rotatedX);
		T translatedX[3];
		translatedX[0] = rotatedX[0]+t[0];
		translatedX[1] = rotatedX[1]+t[1];
		translatedX[2] = rotatedX[2]+t[2];
		T x,y;
		x = k[0]*translatedX[0]+k[4]*translatedX[1]+k[1]*translatedX[2];
		y = k[3]*k[0]*translatedX[1]+k[2]*translatedX[2];
		residuals[0] = T(weight)*(T(observed_x) - x / translatedX[2]);
		residuals[1] = T(weight)*(T(observed_y) - y / translatedX[2]);
		return true;
	}

	private:
	double observed_x;
	double observed_y;
	double k[5];
    double q0[4];
	double weight;
};

// 3 - skew and aspect ratio is fixed
// 9 parameters
// block - [fx,x0,y0,q1,q2,q3,t0,t1,t2]
class ReprojectionErrorFxX0Y0QT{
	public:
		ReprojectionErrorFxX0Y0QT(double observed_x, double observed_y, double * q0, double s, double ar, double weight)
	:observed_x(observed_x),observed_y(observed_y){
		this->q0[0] = q0[0];
		this->q0[1] = q0[1];
		this->q0[2] = q0[2];
		this->q0[3] = q0[3];
		this->ar = ar;
		this->s = s;
		this->weight = weight;
	}
	

	template <typename T>
	bool operator()(const T* const kqt, const T* const X, T* residuals)const{
		const T * k = kqt+6;
		const T * q = kqt;
		const T * t = kqt+3;
		T w = ceres::sqrt(T(1)-q[0]*q[0]-q[1]*q[1]-q[2]*q[2]);
		const T qlocal[4]= {w,q[0],q[1],q[2]};
		T qtotal[4];
		const T tq0[4] = {T(q0[0]), T(q0[1]), T(q0[2]),T(q0[3])};
		QuaternionProduct( qlocal, tq0, qtotal);
		T rotatedX[3];
		QuaternionRotatePoint(qtotal, X, rotatedX);
		T translatedX[3];
		translatedX[0] = rotatedX[0]+t[0];
		translatedX[1] = rotatedX[1]+t[1];
		translatedX[2] = rotatedX[2]+t[2];
		T x,y;
		x = k[0]*translatedX[0]+s*translatedX[1]+k[1]*translatedX[2];
		y = ar*k[0]*translatedX[1]+k[2]*translatedX[2];
		residuals[0] = T(weight)*(T(observed_x) - x / translatedX[2]);
		residuals[1] = T(weight)*(T(observed_y) - y / translatedX[2]);
		return true;
	}

	private:
	double observed_x;
	double observed_y;
	double s,ar;
    double q0[4];
	double weight;
};

// 4 - skew is fixed
// 10 parameters
// block - [fx,x0,y0,ar,q1,q2,q3,t0,t1,t2]
class ReprojectionErrorArFxX0Y0QT{
	public:
		ReprojectionErrorArFxX0Y0QT(double observed_x, double observed_y, double * q0, double s, double weight)
	:observed_x(observed_x),observed_y(observed_y){
		this->q0[0] = q0[0];
		this->q0[1] = q0[1];
		this->q0[2] = q0[2];
		this->q0[3] = q0[3];
		this->s = s;
		this->weight = weight;
	}
	

	template <typename T>
	bool operator()(const T* const kqt, const T* const X, T* residuals)const{
		const T * k = kqt+6;
		const T * q = kqt;
		const T * t = kqt+3;
		T w = ceres::sqrt(T(1)-q[0]*q[0]-q[1]*q[1]-q[2]*q[2]);
		const T qlocal[4]= {w,q[0],q[1],q[2]};
		T qtotal[4];
		const T tq0[4] = {T(q0[0]), T(q0[1]), T(q0[2]),T(q0[3])};
		QuaternionProduct( qlocal, tq0, qtotal);
		T rotatedX[3];
		QuaternionRotatePoint(qtotal, X, rotatedX);
		T translatedX[3];
		translatedX[0] = rotatedX[0]+t[0];
		translatedX[1] = rotatedX[1]+t[1];
		translatedX[2] = rotatedX[2]+t[2];
		T x,y;
		x = k[0]*translatedX[0]+s*translatedX[1]+k[1]*translatedX[2];
		y = k[3]*k[0]*translatedX[1]+k[2]*translatedX[2];
		residuals[0] = T(weight)*(T(observed_x) - x / translatedX[2]);
		residuals[1] = T(weight)*(T(observed_y) - y / translatedX[2]);
		return true;
	}

	private:
	double observed_x;
	double observed_y;
	double s;
    double q0[4];
	double weight;
};

// 5 - all calibration except focal length is fixed,
// 7 parameters
// block - [f,q1,q2,q3,t0,t1,t2]
class ReprojectionErrorFQT{
	public:
		ReprojectionErrorFQT(double observed_x, double observed_y, double * q0, double *k, double weight)
	:observed_x(observed_x),observed_y(observed_y){
		this->q0[0] = q0[0];
		this->q0[1] = q0[1];
		this->q0[2] = q0[2];
		this->q0[3] = q0[3];
		this->k[0] = k[0];
		this->k[1] = k[1];
		this->k[2] = k[2];
		this->k[3] = k[3];
		this->weight = weight;
	}
	

	template <typename T>
	bool operator()(const T* const kqt, const T* const X, T* residuals)const{
		const T f = kqt[6];
		const T * q = kqt;
		const T * t = kqt+3;
		T w = ceres::sqrt(T(1)-q[0]*q[0]-q[1]*q[1]-q[2]*q[2]);
		const T qlocal[4]= {w,q[0],q[1],q[2]};
		T qtotal[4];
		const T tq0[4] = {T(q0[0]), T(q0[1]), T(q0[2]),T(q0[3])};
		QuaternionProduct( qlocal, tq0, qtotal);
		T rotatedX[3];
		QuaternionRotatePoint(qtotal, X, rotatedX);
		T translatedX[3];
		translatedX[0] = rotatedX[0]+t[0];
		translatedX[1] = rotatedX[1]+t[1];
		translatedX[2] = rotatedX[2]+t[2];
		T x,y;
		x = f*translatedX[0]+k[3]*translatedX[1]+k[0]*translatedX[2];
		y = k[2]*f*translatedX[1]+k[1]*translatedX[2];
		residuals[0] = T(weight)*(T(observed_x) - x / translatedX[2]);
		residuals[1] = T(weight)*(T(observed_y) - y / translatedX[2]);
		return true;
	}

	private:
	double observed_x;
	double observed_y;
	double s;
	double k[4];
    double q0[4];
	double weight;
};

//versions with radial distortion

// 1 - all parameters adjusted
// 13 parameters
// block - [fx,x0,y0,ar,s,r0,r1,q1,q2,q3,t0,t1,t2]
class ReprojectionErrorKQTRadial{
	public:
		ReprojectionErrorKQTRadial(double observed_x, double observed_y, double * q0, double weight)
	:observed_x(observed_x),observed_y(observed_y){
		this->q0[0] = q0[0];
		this->q0[1] = q0[1];
		this->q0[2] = q0[2];
		this->q0[3] = q0[3];
		this->weight = weight;
	}
	

	template <typename T>
	bool operator()(const T* const kqt, const T* const X, T* residuals)const{
		const T * k = kqt+8;
		const T * r = kqt+6;
		const T * q = kqt;
		const T * t = kqt+3;
		T w = ceres::sqrt(T(1)-q[0]*q[0]-q[1]*q[1]-q[2]*q[2]);
		const T qlocal[4]= {w,q[0],q[1],q[2]};
		T qtotal[4];
		const T tq0[4] = {T(q0[0]), T(q0[1]), T(q0[2]),T(q0[3])};
		QuaternionProduct( qlocal, tq0, qtotal);
		T rotatedX[3];
		QuaternionRotatePoint(qtotal, X, rotatedX);
		T translatedX[3];
		translatedX[0] = rotatedX[0]+t[0];
		translatedX[1] = rotatedX[1]+t[1];
		translatedX[2] = rotatedX[2]+t[2];
		T x,y;
		x = translatedX[0]/translatedX[2];
		y = translatedX[1]/translatedX[2];
		T p_norm2 = x*x+y*y;
		T dist = T(1)+r[0]*p_norm2+r[1]*p_norm2*p_norm2;
		x=dist*x;
		y=dist*y;
		x = k[0]*x+k[4]*y+k[1];
		y = k[3]*k[0]*y+k[2];
		residuals[0] = T(weight)*(T(observed_x) - x);
		residuals[1] = T(weight)*(T(observed_y) - y);
		return true;
	}

	private:
	double observed_x;
	double observed_y;
    double q0[4];
	double weight;
};

// 2 - K is fixed, provided calibration is used and kept constant
// 6 parameters
// block - [q1,q2,q3,t0,t1,t2]
class ReprojectionErrorQTRadial{
	public:
		ReprojectionErrorQTRadial(double observed_x, double observed_y, double * q0, double * k, double *r, double weight)
	:observed_x(observed_x),observed_y(observed_y){
		this->q0[0] = q0[0];
		this->q0[1] = q0[1];
		this->q0[2] = q0[2];
		this->q0[3] = q0[3];
		this->k[0] = k[0];
		this->k[1] = k[1];
		this->k[2] = k[2];
		this->k[3] = k[3];
		this->k[4] = k[4];
		this->r[0] = r[0];
		this->r[1] = r[1];
		this->weight = weight;
	}
	

	template <typename T>
	bool operator()(const T* const qt, const T* const X, T* residuals)const{
		const T * q = qt;
		const T * t = qt+3;
		T w = ceres::sqrt(T(1)-q[0]*q[0]-q[1]*q[1]-q[2]*q[2]);
		const T qlocal[4]= {w,q[0],q[1],q[2]};
		T qtotal[4];
		const T tq0[4] = {T(q0[0]), T(q0[1]), T(q0[2]),T(q0[3])};
		QuaternionProduct( qlocal, tq0, qtotal);
		T rotatedX[3];
		QuaternionRotatePoint(qtotal, X, rotatedX);
		T translatedX[3];
		translatedX[0] = rotatedX[0]+t[0];
		translatedX[1] = rotatedX[1]+t[1];
		translatedX[2] = rotatedX[2]+t[2];
		T x,y;
		x = translatedX[0]/translatedX[2];
		y = translatedX[1]/translatedX[2];
		T p_norm2 = x*x+y*y;
		T dist = T(1)+r[0]*p_norm2+r[1]*p_norm2*p_norm2;
		x=dist*x;
		y=dist*y;
		x = k[0]*x+k[4]*y+k[1];
		y = k[3]*k[0]*y+k[2];
		residuals[0] = T(weight)*(T(observed_x) - x);
		residuals[1] = T(weight)*(T(observed_y) - y);
		return true;
	}

	private:
	double observed_x;
	double observed_y;
	double k[5],r[2];
    double q0[4];
	double weight;
};

// 3 - skew and aspect ratio is fixed
// 11 parameters
// block - [fx,x0,y0,r0,r1,q1,q2,q3,t0,t1,t2]
class ReprojectionErrorFxX0Y0QTRadial{
	public:
		ReprojectionErrorFxX0Y0QTRadial(double observed_x, double observed_y, double * q0, double s, double ar, double weight)
	:observed_x(observed_x),observed_y(observed_y){
		this->q0[0] = q0[0];
		this->q0[1] = q0[1];
		this->q0[2] = q0[2];
		this->q0[3] = q0[3];
		this->ar = ar;
		this->s = s;
		this->weight = weight;
	}
	

	template <typename T>
	bool operator()(const T* const kqt, const T* const X, T* residuals)const{
		const T *  k = kqt+8;
		const T *  r = kqt+6;
		const T *  q = kqt;
		const T *  t = kqt+3;
		T w = ceres::sqrt(T(1)-q[0]*q[0]-q[1]*q[1]-q[2]*q[2]);
		const T qlocal[4]= {w,q[0],q[1],q[2]};
		T qtotal[4];
		const T tq0[4] = {T(q0[0]), T(q0[1]), T(q0[2]),T(q0[3])};
		QuaternionProduct( qlocal, tq0, qtotal);
		T rotatedX[3];
		QuaternionRotatePoint(qtotal, X, rotatedX);
		T translatedX[3];
		translatedX[0] = rotatedX[0]+t[0];
		translatedX[1] = rotatedX[1]+t[1];
		translatedX[2] = rotatedX[2]+t[2];
		T x,y;
		x = translatedX[0]/translatedX[2];
		y = translatedX[1]/translatedX[2];
		T p_norm2 = x*x+y*y;
		T dist = T(1)+r[0]*p_norm2+r[1]*p_norm2*p_norm2;
		x=dist*x;
		y=dist*y;
		x = k[0]*x+s*y+k[1];
		y = ar*k[0]*y+k[2];
		residuals[0] = T(weight)*(T(observed_x) - x);
		residuals[1] = T(weight)*(T(observed_y) - y);
		return true;
	}

	private:
	double observed_x;
	double observed_y;
	double s,ar;
    double q0[4];
	double weight;
};

// 4 - skew is fixed
// 12 parameters
// block - [fx,x0,y0,ar,r0,r1,q1,q2,q3,t0,t1,t2]
class ReprojectionErrorArFxX0Y0QTRadial{
	public:
		ReprojectionErrorArFxX0Y0QTRadial(double observed_x, double observed_y, double * q0, double s, double weight)
	:observed_x(observed_x),observed_y(observed_y){
		this->q0[0] = q0[0];
		this->q0[1] = q0[1];
		this->q0[2] = q0[2];
		this->q0[3] = q0[3];
		this->s = s;
		this->weight = weight;
	}
	

	template <typename T>
	bool operator()(const T* const kqt, const T* const X, T* residuals)const{
		const T *  k = kqt+8;
		const T *  r = kqt+6;
		const T *  q = kqt;
		const T *  t = kqt+3;
		T w = ceres::sqrt(T(1)-q[0]*q[0]-q[1]*q[1]-q[2]*q[2]);
		const T qlocal[4]= {w,q[0],q[1],q[2]};
		T qtotal[4];
		const T tq0[4] = {T(q0[0]), T(q0[1]), T(q0[2]),T(q0[3])};
		QuaternionProduct( qlocal, tq0, qtotal);
		T rotatedX[3];
		QuaternionRotatePoint(qtotal, X, rotatedX);
		T translatedX[3];
		translatedX[0] = rotatedX[0]+t[0];
		translatedX[1] = rotatedX[1]+t[1];
		translatedX[2] = rotatedX[2]+t[2];
		T x,y;
		x = translatedX[0]/translatedX[2];
		y = translatedX[1]/translatedX[2];
		T p_norm2 = x*x+y*y;
		T dist = T(1)+r[0]*p_norm2+r[1]*p_norm2*p_norm2;
		x=dist*x;
		y=dist*y;
		x = k[0]*x+s*y+k[1];
		y = k[3]*k[0]*y+k[2];
		residuals[0] = T(weight)*(T(observed_x) - x);
		residuals[1] = T(weight)*(T(observed_y) - y);
		return true;
	}

	private:
	double observed_x;
	double observed_y;
	double s;
    double q0[4];
	double weight;
};

// 5 - all calibration except focal length is fixed,
// 9 parameters
// block - [f,r0,r1,q1,q2,q3,t0,t1,t2]
class ReprojectionErrorFQTRadial{
	public:
	ReprojectionErrorFQTRadial(double observed_x, double observed_y, double * q0, double *k, double weight)
	:observed_x(observed_x),observed_y(observed_y){
		this->q0[0] = q0[0];
		this->q0[1] = q0[1];
		this->q0[2] = q0[2];
		this->q0[3] = q0[3];
		this->k[0] = k[0];
		this->k[1] = k[1];
		this->k[2] = k[2];
		this->k[3] = k[3];
		this->weight = weight;
	}
	

	template <typename T>
	bool operator()(const T* const kqt, const T* const X, T* residuals)const{
		const T * fr = kqt+6;
		const T * q = kqt;
		const T * t = kqt+3;
		T w = ceres::sqrt(T(1)-q[0]*q[0]-q[1]*q[1]-q[2]*q[2]);
		const T qlocal[4]= {w,q[0],q[1],q[2]};
		T qtotal[4];
		const T tq0[4] = {T(q0[0]), T(q0[1]), T(q0[2]),T(q0[3])};
		QuaternionProduct( qlocal, tq0, qtotal);
		T rotatedX[3];
		QuaternionRotatePoint(qtotal, X, rotatedX);
		T translatedX[3];
		translatedX[0] = rotatedX[0]+t[0];
		translatedX[1] = rotatedX[1]+t[1];
		translatedX[2] = rotatedX[2]+t[2];
		T x,y;
		x = translatedX[0]/translatedX[2];
		y = translatedX[1]/translatedX[2];
		T p_norm2 = x*x+y*y;
		T dist = T(1)+fr[0]*p_norm2+fr[1]*p_norm2*p_norm2;
		x=dist*x;
		y=dist*y;
		x = fr[2]*x+k[3]*y+k[0];
		y = k[2]*fr[2]*y+k[1];
		residuals[0] = T(weight)*(T(observed_x) - x);
		residuals[1] = T(weight)*(T(observed_y) - y);
		return true;
	}

	private:
	double observed_x;
	double observed_y;
	double k[4];
    double q0[4];
	double weight;
};

//versions with OpenCV distortion model: radial + tangential
// x_corr = x(1+k1*r^2+k2*r^4+k3*r^6)+2*p1*x*y + p2*(r^2+2*x^2)
// y_corr = y(1+k1*r^2+k2*r^4+k3*r^6)+p1*(r^2+2*y^2)+2*p2*x*y
// applied before K

// 3 - skew and aspect ratio is fixed
// 14 parameters
// block - [fx,x0,y0,k1,k2,k3,p1,p2,q1,q2,q3,t0,t1,t2]
class ReprojectionErrorOpenCV{
	public:
		ReprojectionErrorOpenCV(double observed_x, double observed_y, double * q0, double s, double ar, double weight)
	:observed_x(observed_x),observed_y(observed_y){
		this->q0[0] = q0[0];
		this->q0[1] = q0[1];
		this->q0[2] = q0[2];
		this->q0[3] = q0[3];
		this->ar = ar;
		this->s = s;
		this->weight = weight;
	}
	

	template <typename T>
	bool operator()(const T* const kqt, const T* const X, T* residuals)const{
		const T *  k = kqt+11;
		const T *  r = kqt+6;
		const T *  q = kqt;
		const T *  t = kqt+3;
		T w = ceres::sqrt(T(1)-q[0]*q[0]-q[1]*q[1]-q[2]*q[2]);
		const T qlocal[4]= {w,q[0],q[1],q[2]};
		T qtotal[4];
		const T tq0[4] = {T(q0[0]), T(q0[1]), T(q0[2]),T(q0[3])};
		QuaternionProduct( qlocal, tq0, qtotal);
		T rotatedX[3];
		QuaternionRotatePoint(qtotal, X, rotatedX);
		T translatedX[3];
		translatedX[0] = rotatedX[0]+t[0];
		translatedX[1] = rotatedX[1]+t[1];
		translatedX[2] = rotatedX[2]+t[2];
		T x,y;
		x = translatedX[0]/translatedX[2];
		y = translatedX[1]/translatedX[2];
		T p_norm2 = x*x+y*y;
		T radial_dist = T(1)+r[0]*p_norm2+r[1]*p_norm2*p_norm2+r[2]*p_norm2*p_norm2*p_norm2; 
		x=x*radial_dist+T(2)*r[3]*x*y+r[4]*(p_norm2+T(2)*x*x);
		y=y*radial_dist+r[4]*(p_norm2+T(2)*y*x)+T(2)*r[4]*x*y;
		x = k[0]*x+s*y+k[1];
		y = ar*k[0]*y+k[2];
		residuals[0] = T(observed_x) - x;
		residuals[1] = T(observed_y) - y;
		return true;
	}


	private:
	double observed_x;
	double observed_y;
	double s,ar;
    double q0[4];
	double weight;
};




//Projection of calibration boards
class ReprojectionErrorCalib{
	public:
	ReprojectionErrorCalib(double observed_x, double observed_y, double * q0, double s, double ar, double boardx, double boardy,int ptid)
	:observed_x(observed_x),observed_y(observed_y){
		this->q0[0] = q0[0];
		this->q0[1] = q0[1];
		this->q0[2] = q0[2];
		this->q0[3] = q0[3];
		this->ar = ar;
		this->s = s;
		this->boardx = boardx;
		this->boardy = boardy;
		this->ptid=ptid;
	}
	

	template <typename T>
	bool operator()(const T* const kqt, const T* const B, T* residuals)const{
		const T * k = kqt+6;
		const T * q = kqt;
		const T * t = kqt+3;
		const T * eax = B;
		const T * t_board=B+3;
		const T xy_board[3] = {T(boardx),T(boardy),T(0)};
		T rotated_xy_board[3];
		AngleAxisRotatePoint(eax,xy_board,rotated_xy_board);
		T X[3] = {rotated_xy_board[0]+t_board[0],rotated_xy_board[1]+t_board[1],rotated_xy_board[2]+t_board[2]};
		T w = ceres::sqrt(T(1)-q[0]*q[0]-q[1]*q[1]-q[2]*q[2]);
		const T qlocal[4]= {w,q[0],q[1],q[2]};
		T qtotal[4];
		const T tq0[4] = {T(q0[0]), T(q0[1]), T(q0[2]),T(q0[3])};
		QuaternionProduct( qlocal, tq0, qtotal);
		T rotatedX[3];
		QuaternionRotatePoint(qtotal, X, rotatedX);
		T translatedX[3];
		translatedX[0] = rotatedX[0]+t[0];
		translatedX[1] = rotatedX[1]+t[1];
		translatedX[2] = rotatedX[2]+t[2];
		T x,y;
		x = k[0]*translatedX[0]+s*translatedX[1]+k[1]*translatedX[2];
		y = ar*k[0]*translatedX[1]+k[2]*translatedX[2];
		residuals[0] = T(observed_x) - x/translatedX[2];
		residuals[1] = T(observed_y) - y/translatedX[2];
		return true;
	}

	private:
	double observed_x;
	double observed_y;
	double s,ar;
    double q0[4];
	double boardx,boardy;
	int ptid;
};

class ReprojectionErrorMulticamKQT{
	public:
	ReprojectionErrorMulticamKQT(double observed_x, double observed_y, double * qc0,double * qs0)
	:observed_x(observed_x),observed_y(observed_y){
		this->qc0[0] = qc0[0];
		this->qc0[1] = qc0[1];
		this->qc0[2] = qc0[2];
		this->qc0[3] = qc0[3];
		this->qs0[0] = qs0[0];
		this->qs0[1] = qs0[1];
		this->qs0[2] = qs0[2];
		this->qs0[3] = qs0[3];
	}
	

	template <typename T>
	bool operator()(const T* k, const T* qc, const T* tc,const T* qs,const T* ts, const T* const X, T* residuals)const{
		//total relative orientation of the camera
		T w = ceres::sqrt(T(1)-qc[0]*qc[0]-qc[1]*qc[1]-qc[2]*qc[2]);
		const T qclocal[4]= {w,qc[0],qc[1],qc[2]};
		T qctotal[4];
		const T tqc0[4] = {T(qc0[0]), T(qc0[1]), T(qc0[2]),T(qc0[3])};
		QuaternionProduct(tqc0, qclocal, qctotal);
		//total orientation of the system
		w = ceres::sqrt(T(1)-qs[0]*qs[0]-qs[1]*qs[1]-qs[2]*qs[2]);
		const T qslocal[4]= {w,qs[0],qs[1],qs[2]};
		T qstotal[4];
		const T tqs0[4] = {T(qs0[0]), T(qs0[1]), T(qs0[2]),T(qs0[3])};
		QuaternionProduct(tqs0, qslocal, qstotal);
		//total rotation
		T qtotal[4];
		QuaternionProduct(qctotal, qstotal, qtotal);
		//t total = Rc*ts+tc
		T temp[3];
		QuaternionRotatePoint(qctotal,ts,temp);
		T t[3];
		t[0]=temp[0]+tc[0];
		t[1]=temp[1]+tc[1];
		t[2]=temp[2]+tc[2];

		T rotatedX[3];
		QuaternionRotatePoint(qtotal, X, rotatedX);
		T translatedX[3];
		translatedX[0] = rotatedX[0]+t[0];
		translatedX[1] = rotatedX[1]+t[1];
		translatedX[2] = rotatedX[2]+t[2];
		T x,y;
		x = k[0]*translatedX[0]+k[4]*translatedX[1]+k[1]*translatedX[2];
		y = k[3]*k[0]*translatedX[1]+k[2]*translatedX[2];
		residuals[0] = T(observed_x) - x/translatedX[2];
		residuals[1] = T(observed_y) - y/translatedX[2];
		return true;
	}

	private:
	double observed_x;
	double observed_y;
    double qc0[4];
	double qs0[4];
};

// BUNDLER V04 with radial
// 12 parameters
// block - [fx,x0,y0,ar,r0,r1,q1,q2,q3,t0,t1,t2]
class ReprojectionErrorBundlerRadial{
	public:
	ReprojectionErrorBundlerRadial(double observed_x, double observed_y, double * q0, double *k)
	:observed_x(observed_x),observed_y(observed_y){
		this->q0[0] = q0[0];
		this->q0[1] = q0[1];
		this->q0[2] = q0[2];
		this->q0[3] = q0[3];
		this->k[0] = k[0];
		this->k[1] = k[1];
		this->k[2] = k[2];
		this->k[3] = k[3];
	}
	

	template <typename T>
	bool operator()(const T* const kqt, const T* const X, T* residuals)const{
		const T *  fr = kqt+6;
		const T *  r = kqt+7;
		const T *  q = kqt;
		const T *  C = kqt+3;
		T w = ceres::sqrt(T(1)-q[0]*q[0]-q[1]*q[1]-q[2]*q[2]);
		const T qlocal[4]= {w,q[0],q[1],q[2]};
		T qtotal[4];
		const T tq0[4] = {T(q0[0]), T(q0[1]), T(q0[2]),T(q0[3])};
		QuaternionProduct( qlocal, tq0, qtotal);
		T translatedX[3];
		translatedX[0] = X[0]-C[0];
		translatedX[1] = X[1]-C[1];
		translatedX[2] = X[2]-C[2];
		T rotatedX[3];
		QuaternionRotatePoint(qtotal, translatedX, rotatedX);
		T x,y;
		x = rotatedX[0]/rotatedX[2];
		y = rotatedX[1]/rotatedX[2];	
		T p_norm2 = x*x+y*y;
		T dist = T(1)+fr[0]*p_norm2+fr[1]*p_norm2*p_norm2;
		x=dist*x;
		y=dist*y;
		x = fr[2]*x+k[3]*y+k[0];
		y = k[2]*fr[2]*y+k[1];
		residuals[0] = T(observed_x) - x;
		residuals[1] = T(observed_y) - y;
		return true;
	}

	private:
	double observed_x;
	double observed_y;
	double k[4];
    double q0[4];
};

// BUNDLER V04 with radial and shared camera calibration
// 12 parameters
// block - [fx,x0,y0,ar,r0,r1,q1,q2,q3,t0,t1,t2]
class ReprojectionErrorBundlerRadialShared{
	public:
	ReprojectionErrorBundlerRadialShared(double observed_x, double observed_y, double * q0, double *k)
	:observed_x(observed_x),observed_y(observed_y){
		this->q0[0] = q0[0];
		this->q0[1] = q0[1];
		this->q0[2] = q0[2];
		this->q0[3] = q0[3];
		this->k[0] = k[0];
		this->k[1] = k[1];
		this->k[2] = k[2];
		this->k[3] = k[3];
	}
	

	template <typename T>
	bool operator()(const T* const qt, const T* const K, const T* const X, T* residuals)const{
		const T *  fr = K;
		const T *  r = K+1;
		const T *  q = qt;
		const T *  C = qt+3;
		T w = ceres::sqrt(T(1)-q[0]*q[0]-q[1]*q[1]-q[2]*q[2]);
		const T qlocal[4]= {w,q[0],q[1],q[2]};
		T qtotal[4];
		const T tq0[4] = {T(q0[0]), T(q0[1]), T(q0[2]),T(q0[3])};
		QuaternionProduct( qlocal, tq0, qtotal);
		T translatedX[3];
		translatedX[0] = X[0]-C[0];
		translatedX[1] = X[1]-C[1];
		translatedX[2] = X[2]-C[2];
		T rotatedX[3];
		QuaternionRotatePoint(qtotal, translatedX, rotatedX);
		T x,y;
		x = rotatedX[0]/rotatedX[2];
		y = rotatedX[1]/rotatedX[2];	
		T p_norm2 = x*x+y*y;
		T dist = T(1)+fr[0]*p_norm2+fr[1]*p_norm2*p_norm2;
		x=dist*x;
		y=dist*y;
		x = fr[2]*x+k[3]*y+k[0];
		y = k[2]*fr[2]*y+k[1];
		residuals[0] = T(observed_x) - x;
		residuals[1] = T(observed_y) - y;
		return true;
	}

	private:
	double observed_x;
	double observed_y;
	double k[4];
    double q0[4];
};

// VISUAL SFM with radial and shared camera calibration
// 12 parameters
// block - [fx,x0,y0,ar,r0,r1,q1,q2,q3,t0,t1,t2]
class ReprojectionErrorCalibVisualSfMRadial{
public:
	ReprojectionErrorCalibVisualSfMRadial(double observed_x, double observed_y, double * q0)
		:observed_x(observed_x), observed_y(observed_y){
		this->q0[0] = q0[0];
		this->q0[1] = q0[1];
		this->q0[2] = q0[2];
		this->q0[3] = q0[3];
	}


	template <typename T>
	bool operator()(const T* const kqt, const T* const X, T* residuals)const{
		const T *  k = kqt + 6;
		const T *  q = kqt;
		const T *  C = kqt + 3;
		T w = ceres::sqrt(T(1) - q[0] * q[0] - q[1] * q[1] - q[2] * q[2]);
		const T qlocal[4] = { w, q[0], q[1], q[2] };
		T qtotal[4];
		const T tq0[4] = { T(q0[0]), T(q0[1]), T(q0[2]), T(q0[3]) };
		QuaternionProduct( qlocal, tq0, qtotal);
		T translatedX[3];
		translatedX[0] = X[0] - C[0];
		translatedX[1] = X[1] - C[1];
		translatedX[2] = X[2] - C[2];
		T rotatedX[3];
		QuaternionRotatePoint(qtotal, translatedX, rotatedX);
		T x, y;
		x = rotatedX[0] / rotatedX[2];
		y = rotatedX[1] / rotatedX[2];
		T p_norm2 = x*x + y*y;
		T dist = k[0] * p_norm2;
		x = dist*x;
		y = dist*y;
		x = k[2] * x + k[6] * y + k[3];
		y = k[5] * k[2] * y + k[4];
		residuals[0] = T(observed_x)*dist + x;
		residuals[1] = T(observed_y)*dist + y;
		return true;
	}

private:
	double observed_x;
	double observed_y;
	double q0[4];
};

// VISUAL SFM with radial and shared camera calibration
// 12 parameters
// block - [fx,x0,y0,ar,r0,r1,q1,q2,q3,t0,t1,t2]
class ReprojectionErrorKannala{
public:
	ReprojectionErrorKannala(double observed_x, double observed_y)
		:observed_x(observed_x), observed_y(observed_y){
	}


	template <typename T>
	bool operator()(const T* const p, const T* const X, T* residuals)const{
		const T *  eax = p;
		const T *  t = p+3;
		const T *  k = p+6;
		const T mu = p[8];
		const T mv = p[9];
		const T u0 = p[10];
		const T v0 = p[11];

		T rotatedX[3];
		T translatedX[3];
		AngleAxisRotatePoint(eax, X, rotatedX);
		translatedX[0] = rotatedX[0] + t[0];
		translatedX[1] = rotatedX[1] + t[1];
		translatedX[2] = rotatedX[2] + t[2];
		T theta = ceres::acos(translatedX[2] / ceres::sqrt(translatedX[0] * translatedX[0] + translatedX[1] * translatedX[1] + translatedX[2] * translatedX[2]));
		T r = k[0] * theta + k[1] * theta*theta*theta;
		T x, y;
		T rXc = ceres::sqrt(translatedX[0] * translatedX[0] + translatedX[1] * translatedX[1]);
		if (rXc==T(0)){
			rXc = T(1);
		}
		x = r*translatedX[0] / rXc;
		y = r*translatedX[1] / rXc;
		x = mu * x + u0;
		y = mv * y + v0;
		residuals[0] = T(observed_x) + x;
		residuals[1] = T(observed_y) + y;
		return true;
	}

private:
	double observed_x;
	double observed_y;
};



