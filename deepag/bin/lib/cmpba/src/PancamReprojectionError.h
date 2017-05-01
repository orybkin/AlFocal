//an implementation of a reprojection error function for Pancam
//used parameters are k,q,t
//k = fx,x0

#include "ceres/rotation.h"


class PancamReprojectionErrorKQT{
	public:
	PancamReprojectionErrorKQT(double observed_x, double observed_y, double * q0, double * qc0, double ar, double s)
	:observed_x(observed_x),observed_y(observed_y){
		this->q0 = new double[4];
		this->q0[0] = q0[0];
		this->q0[1] = q0[1];
		this->q0[2] = q0[2];
		this->q0[3] = q0[3];
		this->qc0 = new double[4];
		this->qc0[0] = qc0[0];
		this->qc0[1] = qc0[1];
		this->qc0[2] = qc0[2];
		this->qc0[3] = qc0[3];
		this->ar=ar;
		this->s=s;
	}

	template <typename T>
	bool operator()(const T* const k, const T* const qphi, const T* const Cc, const T* const Qc, const T* const rho, const T* const X, T* residuals)const{
		//calculate rotation
		T w = sqrt(T(1)-qphi[0]*qphi[0]-qphi[1]*qphi[1]-qphi[2]*qphi[2]);
		const T qlocal[4]= {w,qphi[0],qphi[1],qphi[2]};
		T qtotal[4];
		T qcircletotal[4];
		T tq0[4] = {T(q0[0]), T(q0[1]), T(q0[2]),T(q0[3])};
		QuaternionProduct(tq0, qlocal, qtotal);
		T rotatedX[3];
		QuaternionRotatePoint(qtotal, X, rotatedX);
		//calculate circle orientation
		w = sqrt(T(1)-Qc[0]*Qc[0]-Qc[1]*Qc[1]-Qc[2]*Qc[2]);
		const T qclocal[4]= {w,Qc[0],Qc[1],Qc[2]};
		const T tqc0[4] = {T(qc0[0]), T(qc0[1]), T(qc0[2]),T(qc0[3])};
		QuaternionProduct(tqc0, qclocal, qcircletotal);
		//calculate the camera center position
		T translatedX[3];
		T camera_position_circle[3];
		T camera_center[3];
		const T qct[4] =  {qcircletotal[0],qcircletotal[1],qcircletotal[2],qcircletotal[3]};
		const T ci[3] = {rho[0]*cos(qphi[3]),rho[0]*sin(qphi[3]), T(0)};
		QuaternionRotatePoint(qct, ci, camera_position_circle);
		camera_center[0]=camera_position_circle[0]+Cc[0];
		camera_center[1]=camera_position_circle[1]+Cc[1];
		camera_center[2]=camera_position_circle[2]+Cc[2];
		const T cc[3] = {camera_center[0],camera_center[1],camera_center[2]};
		T t[3];
		QuaternionRotatePoint(qtotal, cc, t);
		t[0]=-t[0];
		t[1]=-t[1];
		t[2]=-t[2];

		translatedX[0] = rotatedX[0]+t[0];
		translatedX[1] = rotatedX[1]+t[1];
		translatedX[2] = rotatedX[2]+t[2];
		T x,y;
		x = k[0]*translatedX[0]+k[1]*translatedX[2];
		y = k[0]*translatedX[1]+k[2]*translatedX[2];
		residuals[0] = T(observed_x) - x/translatedX[2];
		residuals[1] = T(observed_y) - y/translatedX[2];
		return true;
	}

	private:
	double observed_x;
	double observed_y;
    double * q0;
    double * qc0;
    double ar,s;


};

class PancamReprojectionErrorKQTRadial{
	public:
	PancamReprojectionErrorKQTRadial(double observed_x, double observed_y, double * q0, double * qc0, double ar, double s)
	:observed_x(observed_x),observed_y(observed_y){
		this->q0 = new double[4];
		this->q0[0] = q0[0];
		this->q0[1] = q0[1];
		this->q0[2] = q0[2];
		this->q0[3] = q0[3];
		this->qc0 = new double[4];
		this->qc0[0] = qc0[0];
		this->qc0[1] = qc0[1];
		this->qc0[2] = qc0[2];
		this->qc0[3] = qc0[3];
		this->ar=ar;
		this->s=s;
	}

	template <typename T>
	bool operator()(const T* const k, const T* const qphi, const T* const Cc, const T* const Qc, const T* const rho, const T* const X, T* residuals)const{
		//calculate rotation
		T w = sqrt(T(1)-qphi[0]*qphi[0]-qphi[1]*qphi[1]-qphi[2]*qphi[2]);
		const T qlocal[4]= {w,qphi[0],qphi[1],qphi[2]};
		T qtotal[4];
		T qcircletotal[4];
		T tq0[4] = {T(q0[0]), T(q0[1]), T(q0[2]),T(q0[3])};
		QuaternionProduct(tq0, qlocal, qtotal);
		T rotatedX[3];
		QuaternionRotatePoint(qtotal, X, rotatedX);
		//calculate circle orientation
		w = sqrt(T(1)-Qc[0]*Qc[0]-Qc[1]*Qc[1]-Qc[2]*Qc[2]);
		const T qclocal[4]= {w,Qc[0],Qc[1],Qc[2]};
		const T tqc0[4] = {T(qc0[0]), T(qc0[1]), T(qc0[2]),T(qc0[3])};
		QuaternionProduct(tqc0, qclocal, qcircletotal);
		//calculate the camera center position
		T translatedX[3];
		T camera_position_circle[3];
		T camera_center[3];
		const T qct[4] =  {qcircletotal[0],qcircletotal[1],qcircletotal[2],qcircletotal[3]};
		const T ci[3] = {rho[0]*cos(qphi[3]),rho[0]*sin(qphi[3]), T(0)};
		QuaternionRotatePoint(qct, ci, camera_position_circle);
		camera_center[0]=camera_position_circle[0]+Cc[0];
		camera_center[1]=camera_position_circle[1]+Cc[1];
		camera_center[2]=camera_position_circle[2]+Cc[2];
		const T cc[3] = {camera_center[0],camera_center[1],camera_center[2]};
		T t[3];
		QuaternionRotatePoint(qtotal, cc, t);
		t[0]=-t[0];
		t[1]=-t[1];
		t[2]=-t[2];

		translatedX[0] = rotatedX[0]+t[0];
		translatedX[1] = rotatedX[1]+t[1];
		translatedX[2] = rotatedX[2]+t[2];
		T x,y;
		x = translatedX[0]/translatedX[2];
		y = translatedX[1]/translatedX[2];
		T p_norm2 = x*x+y*y;
		T r = T(1)+k[5]*p_norm2+k[6]*p_norm2*p_norm2;
		x=r*x;
		y=r*y;
		x = k[0]*x+s*y+k[1];
		y = ar*k[0]*y+k[2];
		residuals[0] = T(observed_x) - x;
		residuals[1] = T(observed_y) - y;
		return true;
	}

	private:
	double observed_x;
	double observed_y;
    double * q0;
    double * qc0;
    double ar,s;


};

