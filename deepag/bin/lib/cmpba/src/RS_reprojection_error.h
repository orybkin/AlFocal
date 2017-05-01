// Rolling Shutter reprojection error functions for different RS models
#include "ceres/rotation.h"


template <typename T>
void vectorLERP(const T * v1, const T * v2, T *v3, T t){
	v3[0]=v1[0]*(T(1)-t)+v2[0]*t;
	v3[1]=v1[1]*(T(1)-t)+v2[1]*t;
	v3[2]=v1[2]*(T(1)-t)+v2[2]*t;
}

template <typename T>
void vectorSLERP(T * v1, T * v2, T *v3, T t){
	T dot = v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2]+v1[3]*v2[3];
	if(dot>T(1)){
		v3[0]=v1[0];
		v3[1]=v1[1];
		v3[2]=v1[2];
		v3[3]=v1[3];
	}else{
		T theta = ceres::acos(dot);
		T sinTheta = ceres::sin(theta);
		if (sinTheta < 1e-8){
			v3[0] = v1[0];
			v3[1] = v1[1];
			v3[2] = v1[2];
			v3[3] = v1[3];
		}
		else{
			T ratio1 = ceres::sin((T(1) - t)*theta) / sinTheta;
			T ratio2 = ceres::sin(t*theta) / sinTheta;
			v3[0] = v1[0] * ratio1 + v2[0] * ratio2;
			v3[1] = v1[1] * ratio1 + v2[1] * ratio2;
			v3[2] = v1[2] * ratio1 + v2[2] * ratio2;
			v3[3] = v1[3] * ratio1 + v2[3] * ratio2;
		}
		
	}
}

template <typename T> 
void X_(const T * v, T *M){
	M[0] = 0; M[1] = -v[3]; M[2] = v[2];
	M[3] = v[3]; M[4] = 0; M[5] = -v[1];
	M[6] = -v[2]; M[7] = v[1]; M[8] = 0;	
}


template <typename T>
void matrixMultiply(const T * R1, const T *R2, T * R3){
	R3[0] = R1[0]*R2[0] + R1[1]*R2[3] + R1[2]*R2[6];
	R3[1] = R1[0]*R2[1] + R1[1]*R2[4] + R1[2]*R2[7];
	R3[2] = R1[0]*R2[2] + R1[1]*R2[5] + R1[2]*R2[8];
	R3[3] = R1[3]*R2[0] + R1[4]*R2[3] + R1[5]*R2[6];
	R3[4] = R1[3]*R2[1] + R1[4]*R2[4] + R1[5]*R2[7];
	R3[5] = R1[3]*R2[2] + R1[4]*R2[5] + R1[5]*R2[8];
	R3[6] = R1[6]*R2[0] + R1[7]*R2[3] + R1[8]*R2[6];
	R3[7] = R1[6]*R2[1] + R1[7]*R2[4] + R1[8]*R2[7];
	R3[8] = R1[6]*R2[2] + R1[7]*R2[5] + R1[8]*R2[8];
}

template <typename T>
void c2R(const T * v, T *R){
	R[0] = (v[0] * v[0] - v[1] * v[1] - v[2] * v[2] + T(1)) / (v[0] * v[0] + v[1] * v[1] + v[2] * v[2] + T(1));
	R[1] = -(T(2)* (v[2] - v[0] * v[1])) / (v[0] * v[0] + v[1] * v[1] + v[2] * v[2] + T(1));
	R[2] = (T(2)* (v[1] + v[0] * v[2])) / (v[0] * v[0] + v[1] * v[1] + v[2] * v[2] + T(1));
	R[3] = (T(2)* (v[2] + v[0] * v[1])) / (v[0] * v[0] + v[1] * v[1] + v[2] * v[2] + T(1));
	R[4] = -(v[0] * v[0] - v[1] * v[1] + v[2] * v[2] - T(1)) / (v[0] * v[0] + v[1] * v[1] + v[2] * v[2] + T(1));
	R[5] = -(T(2)* (v[0] - v[1] * v[2])) / (v[0] * v[0] + v[1] * v[1] + v[2] * v[2] + T(1));
	R[6] = -(T(2)* (v[1] - v[0] * v[2])) / (v[0] * v[0] + v[1] * v[1] + v[2] * v[2] + T(1));
	R[7] = (T(2)* (v[0] + v[1] * v[2])) / (v[0] * v[0] + v[1] * v[1] + v[2] * v[2] + T(1));
	R[8] = -(v[0] * v[0] + v[1] * v[1] - v[2] * v[2] - T(1)) / (v[0] * v[0] + v[1] * v[1] + v[2] * v[2] + T(1));
}

template <typename T>
void q2R(const T * q, T *R){
	T s, X, Y, Z, qw, qx, qy, qz,wX, wY, wZ,xX, xY, xZ,yY,yZ,zZ;
	qw = q[0]; qx = q[1]; qy = q[2]; qz = q[3];
	T Nq = qw*qw + qx*qx + qy*qy + qz*qz;	
	CHECK_NE(Nq, T(0));
	 s = T(2) / Nq;
	X = qx*s; Y = qy*s; Z = qz*s;
	wX = qw*X; wY = qw*Y; wZ = qw*Z;
	xX = qx*X; xY = qx*Y; xZ = qx*Z;
	yY = qy*Y; yZ = qy*Z; zZ = qz*Z;
	
	
	R[0] = T(1) - (yY + zZ);
	R[1] = xY - wZ;
	R[2] = xZ + wY;
	R[3] = xY + wZ;
	R[4] = T(1) - (xX + zZ);
	R[5] = yZ - wX;
	R[6] = xZ - wY;
	R[7] = yZ + wX;
	R[8] = T(1) - (xX + yY);
}

// BUNDLER V04 with RS according to Forssen
// 18 parameters (12from first camera and 6 from subsequent)
// block - [fx,x0,y0,ar,r0,r1,q1,q2,q3,t0,t1,t2]
class ReprojectionErrorBundlerRadialRSVideoVertical{
	public:
	ReprojectionErrorBundlerRadialRSVideoVertical(double observed_x, double observed_y, double * q0, double * q02, double *k, double frame_time, double capture_time,double im_height)
	:observed_x(observed_x),observed_y(observed_y){
		this->q0[0] = q0[0];
		this->q0[1] = q0[1];
		this->q0[2] = q0[2];
		this->q0[3] = q0[3];
		this->q02[0] = q02[0];
		this->q02[1] = q02[1];
		this->q02[2] = q02[2];
		this->q02[3] = q02[3];
		this->k[0] = k[0];
		this->k[1] = k[1];
		this->k[2] = k[2];
		this->k[3] = k[3];
		this->capture_time = capture_time;
		this->frame_time = frame_time;
		this->im_height = im_height;
	}
	

	template <typename T>
	bool operator()(const T* const kqt, const T* const kqt2, const T* const X, T* residuals)const{
		const T *  fr = kqt+6;
		const T *  r = kqt+7;
		const T *  q = kqt;
		const T *  C = kqt+3;
		const T * q2 = kqt2;
		const T * C2 = kqt2+3;
		T t = T(im_height)*T(frame_time)/T(observed_y)/T(capture_time);
		T w = ceres::sqrt(T(1)-q[0]*q[0]-q[1]*q[1]-q[2]*q[2]);
		T w2 = ceres::sqrt(T(1)-q2[0]*q2[0]-q2[1]*q2[1]-q2[2]*q2[2]);
		const T qlocal[4]= {w,q[0],q[1],q[2]};
		const T qlocal2[4]= {w2,q2[0],q2[1],q2[2]};
		T qtotal[4];
		T qtotal2[4];
		T qline[4];
		T Cline[3];
		const T tq0[4] = {T(q0[0]), T(q0[1]), T(q0[2]),T(q0[3])};
		const T tq02[4] = {T(q02[0]), T(q02[1]), T(q02[2]),T(q02[3])};
		QuaternionProduct(qlocal, tq0, qtotal);
		QuaternionProduct(qlocal2, tq02, qtotal2);
		vectorSLERP(qtotal,qtotal2,qline,t);
		vectorLERP(C,C2,Cline,t);
		T translatedX[3];
		translatedX[0] = X[0]-Cline[0];
		translatedX[1] = X[1]-Cline[1];
		translatedX[2] = X[2]-Cline[2];
		T rotatedX[3];
		QuaternionRotatePoint(qline, translatedX, rotatedX);
		T x,y;
		x = rotatedX[0]/rotatedX[2];
		y = rotatedX[1]/rotatedX[2];	
		T p_norm2 = x*x+y*y;
		T dist = T(1)+fr[0]*p_norm2+fr[1]*p_norm2*p_norm2;
		x=dist*x;
		y=dist*y;
		x = fr[2]*x+k[3]*y+k[0];
		y = k[2]*fr[2]*y+k[1];
		residuals[0] = T(observed_x) + x;
		residuals[1] = T(observed_y) + y;
		return true;
	}

	private:
	double observed_x;
	double observed_y;
	double k[4];
    double q0[4];
	double q02[4];
	double capture_time;
	double frame_time;
	double im_height;
};

// BUNDLER V04 with RS according to Forssen for horizontal rolling shutter
// 18 parameters (12from first camera and 6 from subsequent)
// block - [fx,x0,y0,ar,r0,r1,q1,q2,q3,t0,t1,t2]
class ReprojectionErrorBundlerRadialRSVideoHorizontal{
	public:
	ReprojectionErrorBundlerRadialRSVideoHorizontal(double observed_x, double observed_y, double * q0, double * q02, double *k, double frame_time, double  capture_time,double im_width)
	:observed_x(observed_x),observed_y(observed_y){
		this->q0[0] = q0[0];
		this->q0[1] = q0[1];
		this->q0[2] = q0[2];
		this->q0[3] = q0[3];
		this->q02[0] = q02[0];
		this->q02[1] = q02[1];
		this->q02[2] = q02[2];
		this->q02[3] = q02[3];
		this->k[0] = k[0];
		this->k[1] = k[1];
		this->k[2] = k[2];
		this->k[3] = k[3];
		this->capture_time = capture_time;
		this->frame_time = frame_time;
		this->im_width = im_width;
	}
	

	template <typename T>
	bool operator()(const T* const kqt, const T* const kqt2, const T* const X, T* residuals)const{
		const T *  fr = kqt+6;
		const T *  r = kqt+7;
		const T *  q = kqt;
		const T *  C = kqt+3;
		const T * q2 = kqt2;
		const T * C2 = kqt2+3;
		//printf("kqt1 = [%f,%f,%f,%f,%f,%f,%f,%f,%f]\n",kqt[0],kqt[1],kqt[2],kqt[3],kqt[4],kqt[5],kqt[6],kqt[7],kqt[8]);
		//printf("kqt2 = [%f,%f,%f,%f,%f,%f,%f,%f,%f]\n",kqt2[0],kqt2[1],kqt2[2],kqt2[3],kqt2[4],kqt2[5],kqt2[6],kqt2[7],kqt2[8]);
		T t = T(observed_y)*T(capture_time)/T(im_width)/T(frame_time);
		T w = ceres::sqrt(T(1)-q[0]*q[0]-q[1]*q[1]-q[2]*q[2]);
		T w2 = ceres::sqrt(T(1)-q2[0]*q2[0]-q2[1]*q2[1]-q2[2]*q2[2]);
		const T qlocal[4]= {w,q[0],q[1],q[2]};
		const T qlocal2[4]= {w2,q2[0],q2[1],q2[2]};
		T qtotal[4];
		T qtotal2[4];
		T qline[4];
		T Cline[3];
		const T tq0[4] = {T(q0[0]), T(q0[1]), T(q0[2]),T(q0[3])};
		const T tq02[4] = {T(q02[0]), T(q02[1]), T(q02[2]),T(q02[3])};
		QuaternionProduct(qlocal, tq0, qtotal);
		QuaternionProduct(qlocal2, tq02, qtotal2);
		vectorSLERP(qtotal,qtotal2,qline,t);
		vectorLERP(C,C2,Cline,t);
		T translatedX[3];
		translatedX[0] = X[0]-Cline[0];
		translatedX[1] = X[1]-Cline[1];
		translatedX[2] = X[2]-Cline[2];
		T rotatedX[3];
		QuaternionRotatePoint(qline, translatedX, rotatedX);
		T x,y;
		x = rotatedX[0]/rotatedX[2];
		y = rotatedX[1]/rotatedX[2];	
		T p_norm2 = x*x+y*y;
		T dist = T(1)+fr[0]*p_norm2+fr[1]*p_norm2*p_norm2;
		x=dist*x;
		y=dist*y;
		x = fr[2]*x+k[3]*y+k[0];
		y = k[2]*fr[2]*y+k[1];
		residuals[0] = T(observed_x) + x;
		residuals[1] = T(observed_y) + y;
		return true;
	}

	private:
	double observed_x;
	double observed_y;
	double k[4];
    double q0[4];
	double q02[4];
	double capture_time;
	double frame_time;
	double im_width;
};

class ReprojectionErrorBundlerRadialRSSingle{
public:
	ReprojectionErrorBundlerRadialRSSingle(const double * obsv, double * q0, int direction){
		this->q0[0] = q0[0];
		this->q0[1] = q0[1];
		this->q0[2] = q0[2];
		this->q0[3] = q0[3];
		this->obsv[0] = obsv[0];
		this->obsv[1] = obsv[1];
		this->direction = direction;
	}


	template <typename T>
	bool operator()(const T* const kqt, const T* const v, const T* const w, const T* const X, T* residuals)const{
		const T *  q = kqt;
		const T *  C = kqt + 3;

		//printf("kqt1 = [%f,%f,%f,%f,%f,%f,%f,%f,%f]\n",kqt[0],kqt[1],kqt[2],kqt[3],kqt[4],kqt[5],kqt[6],kqt[7],kqt[8]);
		//printf("kqt2 = [%f,%f,%f,%f,%f,%f,%f,%f,%f]\n",kqt2[0],kqt2[1],kqt2[2],kqt2[3],kqt2[4],kqt2[5],kqt2[6],kqt2[7],kqt2[8]);
		T t = obsv[direction];
		T ww = ceres::sqrt(T(1) - q[0] * q[0] - q[1] * q[1] - q[2] * q[2]);
		const T qlocal[4] = { ww, q[0], q[1], q[2] };
		T qtotal[4];
		T R[9];
		T Rline[9];
		T Cline[3];
		Cline[0] = C[0] + t*v[0];
		Cline[1] = C[1] + t*v[1];
		Cline[2] = C[2] + t*v[2];
		const T tq0[4] = { T(q0[0]), T(q0[1]), T(q0[2]), T(q0[3]) };
		QuaternionProduct(qlocal, tq0, qtotal);
		QuaternionToRotation(qtotal, R);
		T w1 = t*w[0];
		T w2 = t*w[1];
		T w3 = t*w[2];
		T R1_1 = R[0];
		T R1_2 = R[1];
		T R1_3 = R[2];
		T R2_1 = R[3];
		T R2_2 = R[4];
		T R2_3 = R[5];
		T R3_1 = R[6];
		T R3_2 = R[7];
		T R3_3 = R[8];

		Rline[0] = R1_1 - R2_1*w3 + R3_1*w2; Rline[1] = R1_2 - R2_2*w3 + R3_2*w2; Rline[2] = R1_3 - R2_3*w3 + R3_3*w2;
		Rline[3] = R2_1 + R1_1*w3 - R3_1*w1;  Rline[4] = R2_2 + R1_2*w3 - R3_2*w1;  Rline[5] = R2_3 + R1_3*w3 - R3_3*w1;
		Rline[6] = R3_1 - R1_1*w2 + R2_1*w1; Rline[7] = R3_2 - R1_2*w2 + R2_2*w1; Rline[8] = R3_3 - R1_3*w2 + R2_3*w1;
		T X1[3];
		T XX1[3] = {T(0),T(0),T(0)};
		for (int i = 0; i<3; i++){
			for (int j = 0; j<3; j++){
				XX1[i] += Rline[i * 3 + j] * X[j];
			}
		}

		X1[0] = XX1[0] + Cline[0];
		X1[1] = XX1[1] + Cline[1];
		X1[2] = XX1[2] + Cline[2];

		//QuaternionRotatePoint(qline, translatedX, rotatedX);
		/*T cr[3];
		ceres::CrossProduct(X1, XX1, cr);
		T norm = ceres::sqrt((X1[0] * X1[0] + X1[1] * X1[1] + X1[2] * X1[2])*(XX1[0] * XX1[0] + XX1[1] * XX1[1] + XX1[2] * XX1[2]));
		cr[0] = cr[0] / norm;
		cr[1] = cr[1] / norm;
		cr[2] = cr[2] / norm;*/

		/*residuals[0] = ceres::sqrt(cr[0] * cr[0] + cr[1] * cr[1] + cr[2] * cr[2]);*/
		
		residuals[0] = X1[0]/X1[2]-T(obsv[0]);
		residuals[1] = X1[1] / X1[2] - T(obsv[1]);
	
		return true;
	}

private:
	double obsv[2];
	double q0[4];
	int direction;
};

class ReprojectionErrorBundlerRadialRSSingleShared{
public:
	ReprojectionErrorBundlerRadialRSSingleShared(const double * obsv, double * q0, double *k, int direction){
		this->q0[0] = q0[0];
		this->q0[1] = q0[1];
		this->q0[2] = q0[2];
		this->q0[3] = q0[3];
		this->k[0] = k[0];
		this->k[1] = k[1];
		this->k[2] = k[2];
		this->k[3] = k[3];
		this->obsv[0] = obsv[0];
		this->obsv[1] = obsv[1];
		this->direction = direction;
	}


	template <typename T>
	bool operator()(const T* const kqt, const T* const fr, const T* const v, const T* const w, const T* const X, T* residuals)const{
		const T *  q = kqt;
		const T *  C = kqt + 3;
		
		T t;
		if (direction){
			t = T(obsv[1] / (fr[2] * k[2]) - k[1] / (fr[2] * k[2]));
		}
		else{
			t = T(obsv[0] / fr[2] - k[0] / fr[2]);
		}
		T ww = ceres::sqrt(T(1) - q[0] * q[0] - q[1] * q[1] - q[2] * q[2]);
		const T qlocal[4] = { ww, q[0], q[1], q[2] };
		T qtotal[4];
		T R[9];
		T Rline[9];
		T Cline[3];
		Cline[0] = C[0] + t*v[0];
		Cline[1] = C[1] + t*v[1];
		Cline[2] = C[2] + t*v[2];
		const T tq0[4] = { T(q0[0]), T(q0[1]), T(q0[2]), T(q0[3]) };
		QuaternionProduct(qlocal, tq0, qtotal);
		QuaternionToRotation(qtotal, R);
		T w1 = t*w[0];
		T w2 = t*w[1];
		T w3 = t*w[2];
		T R1_1 = R[0];
		T R1_2 = R[1];
		T R1_3 = R[2];
		T R2_1 = R[3];
		T R2_2 = R[4];
		T R2_3 = R[5];
		T R3_1 = R[6];
		T R3_2 = R[7];
		T R3_3 = R[8];

		Rline[0] = R1_1 - R2_1*w3 + R3_1*w2; Rline[1] = R1_2 - R2_2*w3 + R3_2*w2; Rline[2] = R1_3 - R2_3*w3 + R3_3*w2;
		Rline[3] = R2_1 + R1_1*w3 - R3_1*w1;  Rline[4] = R2_2 + R1_2*w3 - R3_2*w1;  Rline[5] = R2_3 + R1_3*w3 - R3_3*w1;
		Rline[6] = R3_1 - R1_1*w2 + R2_1*w1; Rline[7] = R3_2 - R1_2*w2 + R2_2*w1; Rline[8] = R3_3 - R1_3*w2 + R2_3*w1;
		T X1[3];
		T XX1[3] = { T(0), T(0), T(0) };
		for (int i = 0; i<3; i++){
			for (int j = 0; j<3; j++){
				XX1[i] += Rline[i * 3 + j] * X[j];
			}
		}

		X1[0] = XX1[0] + Cline[0];
		X1[1] = XX1[1] + Cline[1];
		X1[2] = XX1[2] + Cline[2];


		T x, y;
		x = X1[0] / X1[2];
		y = X1[1] / X1[2];
		T p_norm2 = x*x + y*y;
		T dist = T(1) + fr[0] * p_norm2 + fr[1] * p_norm2*p_norm2;
		x = dist*x;
		y = dist*y;
		x = fr[2] * x + k[3] * y + k[0];
		y = k[2] * fr[2] * y + k[1];


		residuals[0] = x - T(obsv[0]);
		residuals[1] = y - T(obsv[1]);
		 
		return true;
	}

private:
	double obsv[2];
	double q0[4];
	double k[4];
	int direction;
};

class ReprojectionErrorBundlerRadialRSSingleSharedWithoutWx{
public:
	ReprojectionErrorBundlerRadialRSSingleSharedWithoutWx(const double * obsv, double * q0, double *k, double wx, int direction){
		this->q0[0] = q0[0];
		this->q0[1] = q0[1];
		this->q0[2] = q0[2];
		this->q0[3] = q0[3];
		this->k[0] = k[0];
		this->k[1] = k[1];
		this->k[2] = k[2];
		this->k[3] = k[3];
		this->obsv[0] = obsv[0];
		this->obsv[1] = obsv[1];
		this->direction = direction;
		this->wx = wx;
	}


	template <typename T>
	bool operator()(const T* const kqt, const T* const fr, const T* const v, const T* const w, const T* const X, T* residuals)const{
		const T *  q = kqt;
		const T *  C = kqt + 3;

		T t;
		if (direction){
			t = T(obsv[1] / (fr[2] * k[2]) - k[1] / (fr[2] * k[2]));
		}
		else{
			t = T(obsv[0] / fr[2] - k[0] / fr[2]);
		}
		T ww = ceres::sqrt(T(1) - q[0] * q[0] - q[1] * q[1] - q[2] * q[2]);
		const T qlocal[4] = { ww, q[0], q[1], q[2] };
		T qtotal[4];
		T R[9];
		T Rline[9];
		T Cline[3];
		Cline[0] = C[0] + t*v[0];
		Cline[1] = C[1] + t*v[1];
		Cline[2] = C[2] + t*v[2];
		const T tq0[4] = { T(q0[0]), T(q0[1]), T(q0[2]), T(q0[3]) };
		QuaternionProduct(qlocal, tq0, qtotal);
		QuaternionToRotation(qtotal, R);
		T w1 = t*wx;
		T w2 = t*w[0];
		T w3 = t*w[1];
		T R1_1 = R[0];
		T R1_2 = R[1];
		T R1_3 = R[2];
		T R2_1 = R[3];
		T R2_2 = R[4];
		T R2_3 = R[5];
		T R3_1 = R[6];
		T R3_2 = R[7];
		T R3_3 = R[8];

		Rline[0] = R1_1 - R2_1*w3 + R3_1*w2; Rline[1] = R1_2 - R2_2*w3 + R3_2*w2; Rline[2] = R1_3 - R2_3*w3 + R3_3*w2;
		Rline[3] = R2_1 + R1_1*w3 - R3_1*w1;  Rline[4] = R2_2 + R1_2*w3 - R3_2*w1;  Rline[5] = R2_3 + R1_3*w3 - R3_3*w1;
		Rline[6] = R3_1 - R1_1*w2 + R2_1*w1; Rline[7] = R3_2 - R1_2*w2 + R2_2*w1; Rline[8] = R3_3 - R1_3*w2 + R2_3*w1;
		T X1[3];
		T XX1[3] = { T(0), T(0), T(0) };
		for (int i = 0; i<3; i++){
			for (int j = 0; j<3; j++){
				XX1[i] += Rline[i * 3 + j] * X[j];
			}
		}

		X1[0] = XX1[0] + Cline[0];
		X1[1] = XX1[1] + Cline[1];
		X1[2] = XX1[2] + Cline[2];


		T x, y;
		x = X1[0] / X1[2];
		y = X1[1] / X1[2];
		T p_norm2 = x*x + y*y;
		T dist = T(1) + fr[0] * p_norm2 + fr[1] * p_norm2*p_norm2;
		x = dist*x;
		y = dist*y;
		x = fr[2] * x + k[3] * y + k[0];
		y = k[2] * fr[2] * y + k[1];


		residuals[0] = x - T(obsv[0]);
		residuals[1] = y - T(obsv[1]);

		return true;
	}

private:
	double obsv[2];
	double q0[4];
	double k[4];
	double wx;
	int direction;
};

class ReprojectionErrorBundlerRadialRSSingleWithCalib{
public:
	ReprojectionErrorBundlerRadialRSSingleWithCalib(const double * obsv, double * q0, double *k, int direction)
		{
		this->q0[0] = q0[0];
		this->q0[1] = q0[1];
		this->q0[2] = q0[2];
		this->q0[3] = q0[3];
		this->k[0] = k[0];
		this->k[1] = k[1];
		this->k[2] = k[2];
		this->k[3] = k[3];
		this->obsv[0] = obsv[0];
		this->obsv[1] = obsv[1];
		this->direction = direction;
	}


	template <typename T>
	bool operator()(const T* const kqt, const T* const v, const T* const w, const T* const X, T* residuals)const{
		const T *  q = kqt;
		const T *  C = kqt + 3;
		const T *  fr = kqt + 6;
		T t;
		if (direction){
			t = T(obsv[1] / (fr[2] * k[2]) - k[1] / (fr[2] * k[2]));
		}
		else{
			t = T(obsv[0] / fr[2] - k[0] / fr[2]);
		}
		T ww = ceres::sqrt(T(1) - q[0] * q[0] - q[1] * q[1] - q[2] * q[2]);
		const T qlocal[4] = { ww, q[0], q[1], q[2] };
		T qtotal[4];
		T R[9];
		T Rline[9];
		T Cline[3];
		Cline[0] = C[0] + t*v[0];
		Cline[1] = C[1] + t*v[1];
		Cline[2] = C[2] + t*v[2];
		const T tq0[4] = { T(q0[0]), T(q0[1]), T(q0[2]), T(q0[3]) };
		QuaternionProduct(qlocal, tq0, qtotal);
		QuaternionToRotation(qtotal, R);
		T w1 = t*w[0];
		T w2 = t*w[1];
		T w3 = t*w[2];
		T R1_1 = R[0];
		T R1_2 = R[1];
		T R1_3 = R[2];
		T R2_1 = R[3];
		T R2_2 = R[4];
		T R2_3 = R[5];
		T R3_1 = R[6];
		T R3_2 = R[7];
		T R3_3 = R[8];

		Rline[0] = R1_1 - R2_1*w3 + R3_1*w2; Rline[1] = R1_2 - R2_2*w3 + R3_2*w2; Rline[2] = R1_3 - R2_3*w3 + R3_3*w2;
		Rline[3] = R2_1 + R1_1*w3 - R3_1*w1;  Rline[4] = R2_2 + R1_2*w3 - R3_2*w1;  Rline[5] = R2_3 + R1_3*w3 - R3_3*w1;
		Rline[6] = R3_1 - R1_1*w2 + R2_1*w1; Rline[7] = R3_2 - R1_2*w2 + R2_2*w1; Rline[8] = R3_3 - R1_3*w2 + R2_3*w1;
		T X1[3];
		T XX1[3] = { T(0), T(0), T(0) };
		for (int i = 0; i<3; i++){
			for (int j = 0; j<3; j++){
				XX1[i] += Rline[i * 3 + j] * X[j];
			}
		}

		X1[0] = XX1[0] + Cline[0];
		X1[1] = XX1[1] + Cline[1];
		X1[2] = XX1[2] + Cline[2];


		T x, y;
		x = X1[0] / X1[2];
		y = X1[1] / X1[2];
		T p_norm2 = x*x + y*y;
		T dist = T(1) + fr[0] * p_norm2 + fr[1] * p_norm2*p_norm2;
		x = dist*x;
		y = dist*y;
		x = fr[2] * x + k[3] * y + k[0];
		y = k[2] * fr[2] * y + k[1];


		residuals[0] = x - T(obsv[0]);
		residuals[1] = y - T(obsv[1]);

		return true;
	}

private:
	double obsv[2];
	double q0[4];
	double k[4];
	int direction;
};

class ReprojectionErrorBundlerRadialRSSingleUpVectorShared{
public:
	ReprojectionErrorBundlerRadialRSSingleUpVectorShared(const double * obsv, double * q0, double *k, int  direction)
		{
		this->q0[0] = q0[0];
		this->q0[1] = q0[1];
		this->q0[2] = q0[2];
		this->q0[3] = q0[3];
		this->k[0] = k[0];
		this->k[1] = k[1];
		this->k[2] = k[2];
		this->k[3] = k[3];
		this->obsv[0] = obsv[0];
		this->obsv[1] = obsv[1];
		this->direction = direction;
	}


	template <typename T>
	bool operator()(const T* const kqt, const T* const fr, const T* const v, const T* const w, const T* const q, const T* const X, T* residuals)const{
		const T *  C = kqt;

		T t;
		if (direction){
			t = T(obsv[1] / (fr[2] * k[2]) - k[1] / (fr[2] * k[2]));
		}
		else{
			t = T(obsv[0] / fr[2] - k[0] / fr[2]);
		}

		T Q[9];
		T qr = q[0];

		
	
		Q[0] = (T(1) - qr * qr) / (T(1) + qr * qr);
		Q[1] = T(0);
		Q[2] = T(2) * qr / (T(1) + qr * qr);
		Q[3] = T(0);
		Q[4] = (T(1) + qr * qr) / (T(1) + qr * qr);
		Q[5] = T(0);
		Q[6] = T(-2) * qr / (T(1) + qr * qr);
		Q[7] = T(0);
		Q[8] = (T(1) - qr * qr) / (T(1) + qr * qr);
		
		T qtotal[4];
		T R[9],Rv[9];
		T Rline[9];
		T Cline[3];
		Cline[0] = C[0] + t*v[0];
		Cline[1] = C[1] + t*v[1];
		Cline[2] = C[2] + t*v[2];
		const T tq0[4] = { T(q0[0]), T(q0[1]), T(q0[2]), T(q0[3]) };
		QuaternionToRotation(tq0, Rv);



		//Rv is actually Rv' due to Matrix class storing in row major
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				R[j * 3 + i] = Rv[i*3] * Q[j*3] + Rv[i*3+1] * Q[j*3+1] + Rv[i*3+2] * Q[j*3+2];
			}
		}

		T w1 = t*w[0];
		T w2 = t*w[1];
		T w3 = t*w[2];
		T R1_1 = R[0];
		T R1_2 = R[3];
		T R1_3 = R[6];
		T R2_1 = R[1];
		T R2_2 = R[4];
		T R2_3 = R[7];
		T R3_1 = R[2];
		T R3_2 = R[5];
		T R3_3 = R[8];

		Rline[0] = R1_1 - R2_1*w3 + R3_1*w2; Rline[1] = R1_2 - R2_2*w3 + R3_2*w2; Rline[2] = R1_3 - R2_3*w3 + R3_3*w2;
		Rline[3] = R2_1 + R1_1*w3 - R3_1*w1;  Rline[4] = R2_2 + R1_2*w3 - R3_2*w1;  Rline[5] = R2_3 + R1_3*w3 - R3_3*w1;
		Rline[6] = R3_1 - R1_1*w2 + R2_1*w1; Rline[7] = R3_2 - R1_2*w2 + R2_2*w1; Rline[8] = R3_3 - R1_3*w2 + R2_3*w1;
		T X1[3];
		T XX1[3] = { T(0), T(0), T(0) };
		for (int i = 0; i<3; i++){
			for (int j = 0; j<3; j++){
				XX1[i] += Rline[i * 3 + j] * X[j];
			}
		}

		X1[0] = XX1[0] + Cline[0];
		X1[1] = XX1[1] + Cline[1];
		X1[2] = XX1[2] + Cline[2];


		T x, y;
		x = X1[0] / X1[2];
		y = X1[1] / X1[2];
		T p_norm2 = x*x + y*y;
		T dist = T(1) + fr[0] * p_norm2 + fr[1] * p_norm2*p_norm2;
		x = dist*x;
		y = dist*y;
		x = fr[2] * x + k[3] * y + k[0];
		y = k[2] * fr[2] * y + k[1];


		residuals[0] = x - T(obsv[0]);
		residuals[1] = y - T(obsv[1]);

		return true;
	}

private:
	double obsv[2];
	double q0[4];
	double k[4];
	int direction;
};

class ReprojectionErrorBundlerRadialRSSingleUpVectorVarShared{
public:
	ReprojectionErrorBundlerRadialRSSingleUpVectorVarShared(const double * obsv, double * q0, double *k, int direction){
		this->q0[0] = q0[0];
		this->q0[1] = q0[1];
		this->q0[2] = q0[2];
		this->q0[3] = q0[3];
		this->k[0] = k[0];
		this->k[1] = k[1];
		this->k[2] = k[2];
		this->k[3] = k[3];
		this->obsv[0] = obsv[0];
		this->obsv[1] = obsv[1];
		this->direction = direction;
	}


	template <typename T>
	bool operator()(const T* const kqt, const T* const fr, const T* const v, const T* const w, const T* const qv, const T* const X, T* residuals)const{
		const T *  C = kqt+3;
		const T *  q = kqt;

		T t;
		if (direction){
			t = T(obsv[1] / (fr[2] * k[2]) - k[1] / (fr[2] * k[2]));
		}
		else{
			t = T(obsv[0] / fr[2] - k[0] / fr[2]);
		}

		T Q[9];
		T qr = qv[0];


		T ww = ceres::sqrt(T(1) - q[0] * q[0] - q[1] * q[1] - q[2] * q[2]);
		const T qlocal[4] = { ww, q[0], q[1], q[2] };

		Q[0] = (T(1) - qr * qr) / (T(1) + qr * qr);
		Q[1] = T(0);
		Q[2] = T(2) * qr / (T(1) + qr * qr);
		Q[3] = T(0);
		Q[4] = (T(1) + qr * qr) / (T(1) + qr * qr);
		Q[5] = T(0);
		Q[6] = T(-2) * qr / (T(1) + qr * qr);
		Q[7] = T(0);
		Q[8] = (T(1) - qr * qr) / (T(1) + qr * qr);

		T qtotal[4];
		T R[9], Rv[9];
		T Rline[9];
		T Cline[3];
		Cline[0] = C[0] + t*v[0];
		Cline[1] = C[1] + t*v[1];
		Cline[2] = C[2] + t*v[2];
		const T tq0[4] = { T(q0[0]), T(q0[1]), T(q0[2]), T(q0[3]) };
		QuaternionProduct(qlocal, tq0, qtotal);
		QuaternionToRotation(qtotal, Rv);



		//Rv is actually Rv' due to Matrix class storing in row major
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				R[j * 3 + i] = Rv[i * 3] * Q[j * 3] + Rv[i * 3 + 1] * Q[j * 3 + 1] + Rv[i * 3 + 2] * Q[j * 3 + 2];
			}
		}

		T w1 = t*w[0];
		T w2 = t*w[1];
		T w3 = t*w[2];
		T R1_1 = R[0];
		T R1_2 = R[3];
		T R1_3 = R[6];
		T R2_1 = R[1];
		T R2_2 = R[4];
		T R2_3 = R[7];
		T R3_1 = R[2];
		T R3_2 = R[5];
		T R3_3 = R[8];

		Rline[0] = R1_1 - R2_1*w3 + R3_1*w2; Rline[1] = R1_2 - R2_2*w3 + R3_2*w2; Rline[2] = R1_3 - R2_3*w3 + R3_3*w2;
		Rline[3] = R2_1 + R1_1*w3 - R3_1*w1;  Rline[4] = R2_2 + R1_2*w3 - R3_2*w1;  Rline[5] = R2_3 + R1_3*w3 - R3_3*w1;
		Rline[6] = R3_1 - R1_1*w2 + R2_1*w1; Rline[7] = R3_2 - R1_2*w2 + R2_2*w1; Rline[8] = R3_3 - R1_3*w2 + R2_3*w1;
		T X1[3];
		T XX1[3] = { T(0), T(0), T(0) };
		for (int i = 0; i<3; i++){
			for (int j = 0; j<3; j++){
				XX1[i] += Rline[i * 3 + j] * X[j];
			}
		}

		X1[0] = XX1[0] + Cline[0];
		X1[1] = XX1[1] + Cline[1];
		X1[2] = XX1[2] + Cline[2];


		T x, y;
		x = X1[0] / X1[2];
		y = X1[1] / X1[2];
		T p_norm2 = x*x + y*y;
		T dist = T(1) + fr[0] * p_norm2 + fr[1] * p_norm2*p_norm2;
		x = dist*x;
		y = dist*y;
		x = fr[2] * x + k[3] * y + k[0];
		y = k[2] * fr[2] * y + k[1];


		residuals[0] = x - T(obsv[0]);
		residuals[1] = y - T(obsv[1]);

		return true;
	}

private:
	double obsv[2];
	double q0[4];
	double k[4];
	int direction;
};

class ReprojectionErrorBundlerRadialRSSingleUpVectorWithCalib{
public:
	ReprojectionErrorBundlerRadialRSSingleUpVectorWithCalib(const double * obsv, double * q0, double *k, int direction)
		{
		this->q0[0] = q0[0];
		this->q0[1] = q0[1];
		this->q0[2] = q0[2];
		this->q0[3] = q0[3];
		this->k[0] = k[0];
		this->k[1] = k[1];
		this->k[2] = k[2];
		this->k[3] = k[3];
		this->obsv[0] = obsv[0];
		this->obsv[1] = obsv[1];
		this->direction = direction;
	}


	template <typename T>
	bool operator()(const T* const kqt, const T* const v, const T* const w, const T* const q, const T* const X, T* residuals)const{
		const T *  C = kqt;
		const T *  fr = kqt + 3;
		T qr = q[0];

		T t;
		if (direction){
			t = T(obsv[1] / (fr[2] * k[2]) - k[1] / (fr[2] * k[2]));
		}
		else{
			t = T(obsv[0] / fr[2] - k[0] / fr[2]);
		}

		T Q[9];

		Q[0] = (T(1) - qr * qr) / (T(1) + qr * qr);
		Q[1] = T(0);
		Q[2] = T(2) * qr / (T(1) + qr * qr);
		Q[3] = T(0);
		Q[4] = (T(1) + qr * qr) / (T(1) + qr * qr);
		Q[5] = T(0);
		Q[6] = T(-2) * qr / (T(1) + qr * qr);
		Q[7] = T(0);
		Q[8] = (T(1) - qr * qr) / (T(1) + qr * qr);


		T qtotal[4];
		T R[9], Rv[9];
		T Rline[9];
		T Cline[3];
		Cline[0] = C[0] + t*v[0];
		Cline[1] = C[1] + t*v[1];
		Cline[2] = C[2] + t*v[2];
		const T tq0[4] = { T(q0[0]), T(q0[1]), T(q0[2]), T(q0[3]) };
		QuaternionToRotation(tq0, Rv);

		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				R[j * 3 + i] = Rv[i * 3] * Q[j * 3] + Rv[i * 3 + 1] * Q[j * 3 + 1] + Rv[i * 3 + 2] * Q[j * 3 + 2];
			}
		}
		T w1 = t*w[0];
		T w2 = t*w[1];
		T w3 = t*w[2];
		T R1_1 = R[0];
		T R1_2 = R[3];
		T R1_3 = R[6];
		T R2_1 = R[1];
		T R2_2 = R[4];
		T R2_3 = R[7];
		T R3_1 = R[2];
		T R3_2 = R[5];
		T R3_3 = R[8];

		Rline[0] = R1_1 - R2_1*w3 + R3_1*w2; Rline[1] = R1_2 - R2_2*w3 + R3_2*w2; Rline[2] = R1_3 - R2_3*w3 + R3_3*w2;
		Rline[3] = R2_1 + R1_1*w3 - R3_1*w1;  Rline[4] = R2_2 + R1_2*w3 - R3_2*w1;  Rline[5] = R2_3 + R1_3*w3 - R3_3*w1;
		Rline[6] = R3_1 - R1_1*w2 + R2_1*w1; Rline[7] = R3_2 - R1_2*w2 + R2_2*w1; Rline[8] = R3_3 - R1_3*w2 + R2_3*w1;
		T X1[3];
		T XX1[3] = { T(0), T(0), T(0) };
		for (int i = 0; i<3; i++){
			for (int j = 0; j<3; j++){
				XX1[i] += Rline[i * 3 + j] * X[j];
			}
		}

		X1[0] = XX1[0] + Cline[0];
		X1[1] = XX1[1] + Cline[1];
		X1[2] = XX1[2] + Cline[2];


		T x, y;
		x = X1[0] / X1[2];
		y = X1[1] / X1[2];
		T p_norm2 = x*x + y*y;
		T dist = T(1) + fr[0] * p_norm2 + fr[1] * p_norm2*p_norm2;
		x = dist*x;
		y = dist*y;
		x = fr[2] * x + k[3] * y + k[0];
		y = k[2] * fr[2] * y + k[1];


		residuals[0] = x - T(obsv[0]);
		residuals[1] = y - T(obsv[1]);

		return true;
	}

private:
	double obsv[2];
	double q0[4];
	double k[4];
	int direction;
};

class ReprojectionErrorBundlerRadialRSSingleCayley{
public:
	ReprojectionErrorBundlerRadialRSSingleCayley(const double * obsv, double * q0, int direction)
		{
		this->q0[0] = q0[0];
		this->q0[1] = q0[1];
		this->q0[2] = q0[2];
		this->q0[3] = q0[3];
		this->obsv[0] = obsv[0];
		this->obsv[1] = obsv[1];
		this->direction = direction;
	}


	template <typename T>
	bool operator()(const T* const kqt, const T* const v, const T* const w, const T* const X, T* residuals)const{
		const T *  q = kqt;
		const T *  C = kqt + 3;
		//printf("kqt1 = [%f,%f,%f,%f,%f,%f,%f,%f,%f]\n",kqt[0],kqt[1],kqt[2],kqt[3],kqt[4],kqt[5],kqt[6],kqt[7],kqt[8]);
		//printf("kqt2 = [%f,%f,%f,%f,%f,%f,%f,%f,%f]\n",kqt2[0],kqt2[1],kqt2[2],kqt2[3],kqt2[4],kqt2[5],kqt2[6],kqt2[7],kqt2[8]);
		T t = T(obsv[direction]);
		T ww = ceres::sqrt(T(1) - q[0] * q[0] - q[1] * q[1] - q[2] * q[2]);
		const T qlocal[4] = { ww, q[0], q[1], q[2] };
		T qtotal[4];
		T R[9];
		T Rw[9];
		T Rline[9];
		T Cline[3];
		Cline[0] = C[0] + t*v[0];
		Cline[1] = C[1] + t*v[1];
		Cline[2] = C[2] + t*v[2];
		const T tq0[4] = { T(q0[0]), T(q0[1]), T(q0[2]), T(q0[3]) };
		QuaternionProduct(qlocal, tq0, qtotal);
		QuaternionToRotation(qtotal, R);
		T wt[3];
		wt[0] = t*w[0];
		wt[1] = t*w[1];
		wt[2] = t*w[2];

		c2R(wt, Rw);

		matrixMultiply(Rw, R, Rline);

		T X1[3];
		T XX1[3] = { T(0), T(0), T(0) };
		for (int i = 0; i<3; i++){
			for (int j = 0; j<3; j++){
				XX1[i] += Rline[i * 3 + j] * X[j];
			}
		}

		X1[0] = XX1[0] + Cline[0];
		X1[1] = XX1[1] + Cline[1];
		X1[2] = XX1[2] + Cline[2];

		//QuaternionRotatePoint(qline, translatedX, rotatedX);
		/*T cr[3];
		ceres::CrossProduct(X1, XX1, cr);
		T norm = ceres::sqrt((X1[0] * X1[0] + X1[1] * X1[1] + X1[2] * X1[2])*(XX1[0] * XX1[0] + XX1[1] * XX1[1] + XX1[2] * XX1[2]));
		cr[0] = cr[0] / norm;
		cr[1] = cr[1] / norm;
		cr[2] = cr[2] / norm;*/

		/*residuals[0] = ceres::sqrt(cr[0] * cr[0] + cr[1] * cr[1] + cr[2] * cr[2]);*/
		residuals[0] = X1[0] / X1[2] - T(obsv[0]);
		residuals[1] = X1[1] / X1[2] - T(obsv[1]);
		return true;
	}

private:
	double obsv[2];
	double q0[4];
	int direction;
};

class ReprojectionErrorBundlerRadialRSSingleCayleyShared{
public:
	ReprojectionErrorBundlerRadialRSSingleCayleyShared(const double *obsv, double * q0, double *k, int direction)
		{
		this->q0[0] = q0[0];
		this->q0[1] = q0[1];
		this->q0[2] = q0[2];
		this->q0[3] = q0[3];
		this->k[0] = k[0];
		this->k[1] = k[1];
		this->k[2] = k[2];
		this->k[3] = k[3];
		this->obsv[0] = obsv[0];
		this->obsv[1] = obsv[1];
		this->direction = direction;
	}


	template <typename T>
	bool operator()(const T* const kqt, const T* const fr, const T* const v, const T* const w, const T* const X, T* residuals)const{
		const T *  q = kqt;
		const T *  C = kqt + 3;
		//printf("kqt1 = [%f,%f,%f,%f,%f,%f,%f,%f,%f]\n",kqt[0],kqt[1],kqt[2],kqt[3],kqt[4],kqt[5],kqt[6],kqt[7],kqt[8]);
		//printf("kqt2 = [%f,%f,%f,%f,%f,%f,%f,%f,%f]\n",kqt2[0],kqt2[1],kqt2[2],kqt2[3],kqt2[4],kqt2[5],kqt2[6],kqt2[7],kqt2[8]);
		T t;
		if (direction){
			t = T(obsv[1] / (fr[2] * k[2]) - k[1] / (fr[2] * k[2]));
		}
		else{
			t = T(obsv[0] / fr[2] - k[0] / fr[2]);
		}
		T ww = ceres::sqrt(T(1) - q[0] * q[0] - q[1] * q[1] - q[2] * q[2]);
		const T qlocal[4] = { ww, q[0], q[1], q[2] };
		T qtotal[4];
		T R[9];
		T Rw[9];
		T Rline[9];
		T Cline[3];
		Cline[0] = C[0] + t*v[0];
		Cline[1] = C[1] + t*v[1];
		Cline[2] = C[2] + t*v[2];
		const T tq0[4] = { T(q0[0]), T(q0[1]), T(q0[2]), T(q0[3]) };
		QuaternionProduct(qlocal, tq0, qtotal);
		QuaternionToRotation(qtotal, R);
		T wt[3];
		wt[0] = t*w[0];
		wt[1] = t*w[1];
		wt[2] = t*w[2];

		c2R(wt, Rw);

		matrixMultiply(Rw, R, Rline);
		

		T X1[3];
		T XX1[3] = { T(0), T(0), T(0) };
		for (int i = 0; i<3; i++){
			for (int j = 0; j<3; j++){
				XX1[i] += Rline[i * 3 + j] * X[j];
			}
		}

		X1[0] = XX1[0] + Cline[0];
		X1[1] = XX1[1] + Cline[1];
		X1[2] = XX1[2] + Cline[2];

		T x, y;
		x = X1[0] / X1[2];
		y = X1[1] / X1[2];
		T p_norm2 = x*x + y*y;
		T dist = T(1) + fr[0] * p_norm2 + fr[1] * p_norm2*p_norm2;
		x = dist*x;
		y = dist*y;
		x = fr[2] * x + k[3] * y + k[0];
		y = k[2] * fr[2] * y + k[1];

		//QuaternionRotatePoint(qline, translatedX, rotatedX);
		/*T cr[3];
		ceres::CrossProduct(X1, XX1, cr);
		T norm = ceres::sqrt((X1[0] * X1[0] + X1[1] * X1[1] + X1[2] * X1[2])*(XX1[0] * XX1[0] + XX1[1] * XX1[1] + XX1[2] * XX1[2]));
		cr[0] = cr[0] / norm;
		cr[1] = cr[1] / norm;
		cr[2] = cr[2] / norm;*/

		/*residuals[0] = ceres::sqrt(cr[0] * cr[0] + cr[1] * cr[1] + cr[2] * cr[2]);*/
		residuals[0] = x - T(obsv[0]);
		residuals[1] = y - T(obsv[1]);
		return true;
	}

private:
	double obsv[2];
	double q0[4],k[4];
	int direction;
};

class ReprojectionErrorBundlerRadialRSSingleCayleyWithCalib{
public:
	ReprojectionErrorBundlerRadialRSSingleCayleyWithCalib(const double * obsv, double * q0, double *k, int direction)
		{
		this->q0[0] = q0[0];
		this->q0[1] = q0[1];
		this->q0[2] = q0[2];
		this->q0[3] = q0[3];
		this->k[0] = k[0];
		this->k[1] = k[1];
		this->k[2] = k[2];
		this->k[3] = k[3];
		this->obsv[0] = obsv[0];
		this->obsv[1] = obsv[1];
		this->direction = direction;
	}


	template <typename T>
	bool operator()(const T* const kqt, const T* const v, const T* const w, const T* const X, T* residuals)const{
		const T *  q = kqt;
		const T *  C = kqt + 3;
		const T *  fr = kqt + 6;
		
		T t;
		if (direction){
			t = T(obsv[1] / (fr[2] * k[2]) - k[1] / (fr[2] * k[2]));
		}
		else{
			t = T(obsv[0] / fr[2] - k[0] / fr[2]);
		}
		T ww = ceres::sqrt(T(1) - q[0] * q[0] - q[1] * q[1] - q[2] * q[2]);
		const T qlocal[4] = { ww, q[0], q[1], q[2] };
		T qtotal[4];
		T R[9];
		T Rw[9];
		T Rline[9];
		T Cline[3];
		Cline[0] = C[0] + t*v[0];
		Cline[1] = C[1] + t*v[1];
		Cline[2] = C[2] + t*v[2];
		const T tq0[4] = { T(q0[0]), T(q0[1]), T(q0[2]), T(q0[3]) };
		QuaternionProduct(qlocal, tq0, qtotal);
		QuaternionToRotation(qtotal, R);
		T wt[3];
		wt[0] = t*w[0];
		wt[1] = t*w[1];
		wt[2] = t*w[2];

		c2R(wt, Rw);

		matrixMultiply(Rw, R, Rline);

		T X1[3];
		T XX1[3] = { T(0), T(0), T(0) };
		for (int i = 0; i<3; i++){
			for (int j = 0; j<3; j++){
				XX1[i] += Rline[i * 3 + j] * X[j];
			}
		}

		X1[0] = XX1[0] + Cline[0];
		X1[1] = XX1[1] + Cline[1];
		X1[2] = XX1[2] + Cline[2];

		T x, y;
		x = X1[0] / X1[2];
		y = X1[1] / X1[2];
		T p_norm2 = x*x + y*y;
		T dist = T(1) + fr[0] * p_norm2 + fr[1] * p_norm2*p_norm2;
		x = dist*x;
		y = dist*y;
		x = fr[2] * x + k[3] * y + k[0];
		y = k[2] * fr[2] * y + k[1];

		residuals[0] = x - T(obsv[0]);
		residuals[1] = y - T(obsv[1]);
		return true;
	}

private:
	double obsv[2];
	double q0[4], k[4];
	int direction;
};

// Rolling shutter with linearized translation and SLERP
// q and w are quaternions for the initial and final orientation
// C and t are vectors describing initial and final camera center
// block - [fx,x0,y0,ar,r0,r1,q1,q2,q3,t0,t1,t2]
class ReprojectionErrorRSForssen{
public:
	ReprojectionErrorRSForssen(const double * obsv, double * q0, double *k, int direction)
	{
		this->q0[0] = q0[0];
		this->q0[1] = q0[1];
		this->q0[2] = q0[2];
		this->q0[3] = q0[3];
		this->k[0] = k[0];
		this->k[1] = k[1];
		this->k[2] = k[2];
		this->k[3] = k[3];
		this->direction = direction;
		this->obsv[0] = obsv[0];
		this->obsv[1] = obsv[1];
	}


	template <typename T>
	bool operator()(const T* const kqt, const T* const v, const T* const w, const T* const X, T* residuals)const{
		const T *  fr = kqt + 6;
		const T *  r = kqt + 7;
		const T *  q = kqt;
		const T *  C = kqt + 3;
		const T * q2 = w;
		const T * C2 = v;
		T t;
		if (direction){
			t = T(obsv[1] / (fr[2] * k[2]) - k[1] / (fr[2] * k[2]));
		}
		else{
			t = T(obsv[0] / fr[2] - k[0] / fr[2]);
		}
		T qw = ceres::sqrt(T(1) - q[0] * q[0] - q[1] * q[1] - q[2] * q[2]);
		T qw2 = ceres::sqrt(T(1) - q2[0] * q2[0] - q2[1] * q2[1] - q2[2] * q2[2]);
		const T qlocal[4] = { qw, q[0], q[1], q[2] };
		const T qlocal2[4] = { qw2, q2[0], q2[1], q2[2] };
		T qtotal[4];
		T qtotal2[4];
		T qline[4];
		T Cline[3];
		const T tq0[4] = { T(q0[0]), T(q0[1]), T(q0[2]), T(q0[3]) };
		const T tq02[4] = { T(q0[0]), T(q0[1]), T(q0[2]), T(q0[3]) };
		QuaternionProduct(qlocal, tq0, qtotal);
		QuaternionProduct(qlocal2, tq02, qtotal2);
		vectorSLERP(qtotal, qtotal2, qline, t);
		vectorLERP(C, C2, Cline, t);

		T translatedX[3];
		translatedX[0] = X[0] - Cline[0];
		translatedX[1] = X[1] - Cline[1];
		translatedX[2] = X[2] - Cline[2];
		T rotatedX[3];
		QuaternionRotatePoint(qline, translatedX, rotatedX);
		T x, y;
		x = rotatedX[0] / rotatedX[2];
		y = rotatedX[1] / rotatedX[2];
		T p_norm2 = x*x + y*y;
		T dist = T(1) + fr[0] * p_norm2 + fr[1] * p_norm2*p_norm2;
		x = dist*x;
		y = dist*y;
		x = fr[2] * x + k[3] * y + k[0];
		y = k[2] * fr[2] * y + k[1];
		residuals[0] = T(obsv[0]) - x;
		residuals[1] = T(obsv[1]) - y;
		return true;
	}

private:
	double obsv[2];
	double k[4];
	double q0[4];
	int direction;
};


class ReprojectionErrorBundlerRadialRSSingleEAXWithCalib{
public:
	ReprojectionErrorBundlerRadialRSSingleEAXWithCalib(const double * obsv, double * q0, double *k, int direction)
	{
		this->q0[0] = q0[0];
		this->q0[1] = q0[1];
		this->q0[2] = q0[2];
		this->q0[3] = q0[3];
		this->k[0] = k[0];
		this->k[1] = k[1];
		this->k[2] = k[2];
		this->k[3] = k[3];
		this->obsv[0] = obsv[0];
		this->obsv[1] = obsv[1];
		this->direction = direction;
	}


	template <typename T>
	bool operator()(const T* const kqt, const T* const v, const T* const w, const T* const X, T* residuals)const{
		const T *  q = kqt;
		const T *  C = kqt + 3;
		const T *  fr = kqt + 6;
		T t;
		if (direction){
			t = T(obsv[1] / (fr[2] * k[2]) - k[1] / (fr[2] * k[2]));
		}
		else{
			t = T(obsv[0] / fr[2] - k[0] / fr[2]);
		}
		T ww = ceres::sqrt(T(1) - q[0] * q[0] - q[1] * q[1] - q[2] * q[2]);
		const T qlocal[4] = { ww, q[0], q[1], q[2] };
		T qtotal[4];
		T R[9];
		T Rw[9];
		T Rline[9];
		T Cline[3];
		Cline[0] = C[0] + t*v[0];
		Cline[1] = C[1] + t*v[1];
		Cline[2] = C[2] + t*v[2];
		const T tq0[4] = { T(q0[0]), T(q0[1]), T(q0[2]), T(q0[3]) };
		QuaternionProduct(qlocal, tq0, qtotal);
		QuaternionToRotation(qtotal, R);
		T eax[3];
		T X1[3];
		T XX1[3];
		eax[0] = t*w[0];
		eax[1] = t*w[1];
		eax[2] = t*w[2];
		QuaternionRotatePoint(qtotal, X, X1);
		AngleAxisRotatePoint(eax,X1,XX1);
		
		

		X1[0] = XX1[0] + Cline[0];
		X1[1] = XX1[1] + Cline[1];
		X1[2] = XX1[2] + Cline[2];


		T x, y;
		x = X1[0] / X1[2];
		y = X1[1] / X1[2];
		T p_norm2 = x*x + y*y;
		T dist = T(1) + fr[0] * p_norm2 + fr[1] * p_norm2*p_norm2;
		x = dist*x;
		y = dist*y;
		x = fr[2] * x + k[3] * y + k[0];
		y = k[2] * fr[2] * y + k[1];


		residuals[0] = x - T(obsv[0]);
		residuals[1] = y - T(obsv[1]);

		return true;
	}

private:
	double obsv[2];
	double q0[4];
	double k[4];
	int direction;
};
