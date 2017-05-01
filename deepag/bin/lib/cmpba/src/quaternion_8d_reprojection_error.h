#include "ceres/rotation.h"


//#define DEBUG_STUFF

// Reprojection error according to paper [1]
// parameter block is [q0,q1,q2,q3,t1,t2,t3]

class Quaternion8dReprojectionError{
	public:
	Quaternion8dReprojectionError(double observed_x, double observed_y, double ar, double s)
	:observed_x(observed_x),observed_y(observed_y){
		this->ar=ar;
		this->s=s;
	}

	Quaternion8dReprojectionError(double observed_x, double observed_y, double ar, double s,int camid, int pointid)
	:observed_x(observed_x),observed_y(observed_y){
		this->ar=ar;
		this->s=s;
		this->camid=camid;
		this->pointid=pointid;
	}

	Quaternion8dReprojectionError(double observed_x, double observed_y, double ar, double s, double x0, double y0)
	:observed_x(observed_x),observed_y(observed_y){
		this->ar=ar;
		this->s=s;
		this->camid=camid;
		this->pointid=pointid;
		this->x0=x0;
		this->y0=y0;
	}

	template <typename T>
	bool operator()(const T* kqt, const T* const X, T* residuals)const{
		T q0 = kqt[0];
		T q1 = kqt[1];
		T q2 = kqt[2];
		T q3 = kqt[3];
		T rotatedX[3];
		const T * C = kqt+4;
		//QuaternionRotatePoint(qtotal, X, rotatedX);
		T translatedX[3];
		T R[3][3];
		R[0][0] = (q0 * q0) + (q1 * q1) - (q2 * q2) - (q3 * q3);
		R[0][1] = T(-2) * q0 * q3 + T(2) * q1 * q2;
		R[0][2] = T(2) * q0 * q2 + T(2) * q1 * q3;
		R[1][0] = T(2) * q0 * q3 + T(2) * q1 * q2;
		R[1][1] = q0 * q0 - q1 * q1 + q2 * q2 - q3 * q3;
		R[1][2] = T(-2) * q0 * q1 + T(2) * q2 * q3;
		R[2][0] = T(-2) * q0 * q2 + T(2) * q1 * q3;
		R[2][1] = T(2) * q0 * q1 + T(2) * q2 * q3;
		R[2][2] = q0 * q0 - q1 * q1 - q2 * q2 + q3 * q3;

		translatedX[0] = X[0]-C[0];
		translatedX[1] = X[1]-C[1];
		translatedX[2] = X[2]-C[2];

		rotatedX[0]=R[0][0]*translatedX[0]+R[0][1]*translatedX[1]+R[0][2]*translatedX[2];
		rotatedX[1]=R[1][0]*translatedX[0]+R[1][1]*translatedX[1]+R[1][2]*translatedX[2];
		rotatedX[2]=R[2][0]*translatedX[0]+R[2][1]*translatedX[1]+R[2][2]*translatedX[2];

		T k11 = q0*q0+q1*q1+q2*q2+q3*q3;

		T x,y;
		//x = k[0]*translatedX[0]+k[4]*translatedX[1]+k[1]*translatedX[2];
		//x = k11*rotatedX[0]+s*rotatedX[1]+k[0]*rotatedX[2];
		x = k11*rotatedX[0]+s*rotatedX[1]+x0*rotatedX[2];
		//y = k[0]*k[3]*translatedX[1]+k[2]*translatedX[2];
		//y = ar*k11*rotatedX[1]+k[1]*rotatedX[2];
		y = ar*k11*rotatedX[1]+y0*rotatedX[2];
		residuals[0] = T(observed_x) - x/rotatedX[2];
		residuals[1] = T(observed_y) - y/rotatedX[2];
#ifdef DEBUG_STUFF
		//printf("residuals [x y]: [%f %f] \n",residuals[0],residuals[1]);
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
			//printf("[k q t]: [%f %f %f %f %f %f %f %f %f %f %f] \n",k[0],k[1],k[2],kr[0],kr[1],q[0],q[1],q[2],t[0],t[1],t[2]);
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
    double ar;
    double s,x0,y0;
	int pointid,camid;


};

class Quaternion8dReprojectionErrorRadial{
	public:
	Quaternion8dReprojectionErrorRadial(double observed_x, double observed_y, double ar, double s)
	:observed_x(observed_x),observed_y(observed_y){
		this->ar=ar;
		this->s=s;
	}

	Quaternion8dReprojectionErrorRadial(double observed_x, double observed_y, double ar, double s,int camid, int pointid)
	:observed_x(observed_x),observed_y(observed_y){
		this->ar=ar;
		this->s=s;
		this->camid=camid;
		this->pointid=pointid;
	}

	Quaternion8dReprojectionErrorRadial(double observed_x, double observed_y, double ar, double s, double x0, double y0)
	:observed_x(observed_x),observed_y(observed_y){
		this->ar=ar;
		this->s=s;
		this->camid=camid;
		this->pointid=pointid;
		this->x0=x0;
		this->y0=y0;
	}

	template <typename T>
	bool operator()(const T* kqt, const T* const X, T* residuals)const{
		T q0 = kqt[0];
		T q1 = kqt[1];
		T q2 = kqt[2];
		T q3 = kqt[3];
		T rotatedX[3];
		const T * C = kqt+4;
		//QuaternionRotatePoint(qtotal, X, rotatedX);
		T translatedX[3];
		T R[3][3];
		R[0][0] = (q0 * q0) + (q1 * q1) - (q2 * q2) - (q3 * q3);
		R[0][1] = T(-2) * q0 * q3 + T(2) * q1 * q2;
		R[0][2] = T(2) * q0 * q2 + T(2) * q1 * q3;
		R[1][0] = T(2) * q0 * q3 + T(2) * q1 * q2;
		R[1][1] = q0 * q0 - q1 * q1 + q2 * q2 - q3 * q3;
		R[1][2] = T(-2) * q0 * q1 + T(2) * q2 * q3;
		R[2][0] = T(-2) * q0 * q2 + T(2) * q1 * q3;
		R[2][1] = T(2) * q0 * q1 + T(2) * q2 * q3;
		R[2][2] = q0 * q0 - q1 * q1 - q2 * q2 + q3 * q3;

		translatedX[0] = X[0]-C[0];
		translatedX[1] = X[1]-C[1];
		translatedX[2] = X[2]-C[2];

		rotatedX[0]=R[0][0]*translatedX[0]+R[0][1]*translatedX[1]+R[0][2]*translatedX[2];
		rotatedX[1]=R[1][0]*translatedX[0]+R[1][1]*translatedX[1]+R[1][2]*translatedX[2];
		rotatedX[2]=R[2][0]*translatedX[0]+R[2][1]*translatedX[1]+R[2][2]*translatedX[2];

		T k11 = q0*q0+q1*q1+q2*q2+q3*q3;
		T p_norm2;
		T x,y;
		x = rotatedX[0]/rotatedX[2];
		y = rotatedX[1]/rotatedX[2];
		p_norm2 = x*x+y*y;
		T r = T(1)+kqt[7]*p_norm2+kqt[8]*p_norm2*p_norm2;
		x=r*x;
		y=r*y;
		//x = k[0]*translatedX[0]+k[4]*translatedX[1]+k[1]*translatedX[2];
		//x = k11*x+s*y+k[0];
		x = k11*x+s*y+x0;
		//y = k[0]*k[3]*translatedX[1]+k[2]*translatedX[2];
		//y = ar*k11*y+k[1];
		y = k11*ar*y+y0;
		residuals[0] = T(observed_x) - x;
		residuals[1] = T(observed_y) - y;
#ifdef DEBUG_STUFF
		//printf("residuals [x y]: [%f %f] \n",residuals[0],residuals[1]);
		double res[2];
		double proj[2];
		char cislo[100];		
		sprintf(cislo,"%f",residuals[0]);
		res[0]=atof(cislo);
		sprintf(cislo,"%f",residuals[1]);
		res[1]=atof(cislo);

		if(res[0]>100||res[1]>100){
			sprintf(cislo,"%f",x);
			proj[0]=atof(cislo);
			sprintf(cislo,"%f",y);
			proj[1]=atof(cislo);
			//printf("[k q t]: [%f %f %f %f %f %f %f %f %f %f %f] \n",k[0],k[1],k[2],kr[0],kr[1],q[0],q[1],q[2],t[0],t[1],t[2]);
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
    double ar;
    double s,x0,y0;
	int pointid,camid;


};



