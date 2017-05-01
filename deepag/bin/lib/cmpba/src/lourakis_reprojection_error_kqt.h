//an implementation of a reprojection error function according to Lourakis
//used parameters are k,q,t
//k = fx,x0

#include "ceres/rotation.h"



namespace ceres {
namespace examples {

//#define DEBUG_STUFF


class LourakisReprojectionErrorKQT{
	public:
	LourakisReprojectionErrorKQT(double observed_x, double observed_y, double * q0, double ar, double s)
	:observed_x(observed_x),observed_y(observed_y){
		this->q0 = new double[4];
		this->q0[0] = q0[0];
		this->q0[1] = q0[1];
		this->q0[2] = q0[2];
		this->q0[3] = q0[3];
		this->ar=ar;
		this->s=s;
	}

	LourakisReprojectionErrorKQT(double observed_x, double observed_y, double * q0, double ar, double s,int camid, int pointid)
	:observed_x(observed_x),observed_y(observed_y){
		this->q0 = new double[4];
		this->q0[0] = q0[0];
		this->q0[1] = q0[1];
		this->q0[2] = q0[2];
		this->q0[3] = q0[3];
		this->ar=ar;
		this->s=s;
		this->camid=camid;
		this->pointid=pointid;
	}

	LourakisReprojectionErrorKQT(double observed_x, double observed_y, double * q0, double ar, double s, double x0, double y0)
	:observed_x(observed_x),observed_y(observed_y){
		this->q0 = new double[4];
		this->q0[0] = q0[0];
		this->q0[1] = q0[1];
		this->q0[2] = q0[2];
		this->q0[3] = q0[3];
		this->ar=ar;
		this->s=s;
		this->x0=x0;
		this->y0=y0;
		this->camid=camid;
		this->pointid=pointid;
	}

	template <typename T>
	bool operator()(const T* const k, const T* const q, const T* const t, const T* const X, T* residuals)const{
		T w = sqrt(T(1)-q[0]*q[0]-q[1]*q[1]-q[2]*q[2]);
		const T qlocal[4]= {w,q[0],q[1],q[2]};
		T qtotal[4];
		const T tq0[4] = {T(q0[0]), T(q0[1]), T(q0[2]),T(q0[3])};
		QuaternionProduct(tq0, qlocal, qtotal);
		T rotatedX[3];
		QuaternionRotatePoint(qtotal, X, rotatedX);
		T translatedX[3];
		translatedX[0] = rotatedX[0]+t[0];
		translatedX[1] = rotatedX[1]+t[1];
		translatedX[2] = rotatedX[2]+t[2];
		T x,y;
		//x = k[0]*translatedX[0]+k[4]*translatedX[1]+k[1]*translatedX[2];
		//x = k[0]*translatedX[0]+s*translatedX[1]+k[1]*translatedX[2];
		x = k[0]*translatedX[0]+s*translatedX[1]+x0*translatedX[2];
		//y = k[0]*k[3]*translatedX[1]+k[2]*translatedX[2];
		//y = ar*k[0]*translatedX[1]+k[2]*translatedX[2];
		y = ar*k[0]*translatedX[1]+y0*translatedX[2];
		residuals[0] = T(observed_x) - x/translatedX[2];
		residuals[1] = T(observed_y) - y/translatedX[2];
#ifdef DEBUG_STUFF
		printf("residuals [x y]: [%f %f] \n",residuals[0],residuals[1]);
#endif
		return true;
	}

	private:
	double observed_x;
	double observed_y;
    double * q0;
    double ar;
    double s,x0,y0;
	int pointid,camid;


};

class LourakisReprojectionErrorKQTRadial{
	public:
	LourakisReprojectionErrorKQTRadial(double observed_x, double observed_y, double * q0, double ar, double s)
	:observed_x(observed_x),observed_y(observed_y){
		this->q0 = new double[4];
		this->q0[0] = q0[0];
		this->q0[1] = q0[1];
		this->q0[2] = q0[2];
		this->q0[3] = q0[3];
		this->ar=ar;
		this->s=s;
	}

	LourakisReprojectionErrorKQTRadial(double observed_x, double observed_y, double * q0, double ar, double s, int camid, int pointid)
	:observed_x(observed_x),observed_y(observed_y){
		this->q0 = new double[4];
		this->q0[0] = q0[0];
		this->q0[1] = q0[1];
		this->q0[2] = q0[2];
		this->q0[3] = q0[3];
		this->ar=ar;
		this->s=s;
		this->camid=camid;
		this->pointid=pointid;
	}

	LourakisReprojectionErrorKQTRadial(double observed_x, double observed_y, double * q0, double ar, double s, double x0, double y0)
	:observed_x(observed_x),observed_y(observed_y){
		this->q0 = new double[4];
		this->q0[0] = q0[0];
		this->q0[1] = q0[1];
		this->q0[2] = q0[2];
		this->q0[3] = q0[3];
		this->ar=ar;
		this->s=s;
		this->x0=x0;
		this->y0=y0;
		this->camid=camid;
		this->pointid=pointid;
	}

	template <typename T>
	bool operator()(const T* const k,const T* const kr, const T* const q, const T* const t, const T* const X, T* residuals)const{
		T w = sqrt(T(1)-q[0]*q[0]-q[1]*q[1]-q[2]*q[2]);
		const T qlocal[4]= {w,q[0],q[1],q[2]};
		T qtotal[4];
		const T tq0[4] = {T(q0[0]), T(q0[1]), T(q0[2]),T(q0[3])};
		QuaternionProduct(tq0, qlocal, qtotal);
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
		T r = T(1)+kr[0]*p_norm2+kr[1]*p_norm2*p_norm2;
		x=r*x;
		y=r*y;
		//x = k[0]*translatedX[0]+k[4]*translatedX[1]+k[1]*translatedX[2];
		//x = k[0]*x+s*y+k[1];
		x = k[0]*x+s*y+x0;
		//y = k[0]*k[3]*translatedX[1]+k[2]*translatedX[2];
		//y = ar*k[0]*y+k[2];
		y = ar*k[0]*y+y0;
		residuals[0] = T(observed_x) - x;
		residuals[1] = T(observed_y) - y;
#ifdef DEBUG_STUFF
		double res[2];
		double proj[2];
		char cislo[100];		
		sprintf(cislo,"%f",residuals[0]);
		res[0]=atof(cislo);
		sprintf(cislo,"%f",residuals[1]);
		res[1]=atof(cislo);

		if(res[0]>1000||res[1]>1000){
			sprintf(cislo,"%f",x);
			proj[0]=atof(cislo);
			sprintf(cislo,"%f",y);
			proj[1]=atof(cislo);
			printf("[k q t]: [%f %f %f %f %f %f %f %f %f %f %f] \n",k[0],k[1],k[2],kr[0],kr[1],q[0],q[1],q[2],t[0],t[1],t[2]);
			printf("[X]: [%f %f %f ] \n",X[0],X[1],X[2]);
			printf("residuals [x y]: [%f %f] \n",residuals[0],residuals[1]);
			printf("observations [x y]: [%f %f] \n",observed_x,observed_y);
			printf("projection [x y]: [%f %f] \n",proj[0],proj[1]);
			printf("camera id = %d, point id =%d",camid,pointid);
		}
#endif
		return true;
	}

	private:
	double observed_x;
	double observed_y;
    double * q0;
    double ar;
    double s,x0,y0;
	int pointid,camid;


};

}}
