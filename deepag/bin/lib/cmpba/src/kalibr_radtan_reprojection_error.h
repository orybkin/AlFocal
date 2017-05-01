#ifndef CERES_EXAMPLES_RADTAN_REPROJECTION_ERROR_H_
#define CERES_EXAMPLES_RADTAN_REPROJECTION_ERROR_H_

#include "ceres/rotation.h"

using namespace ceres;



struct KalibrRadtanReprojectionError {
	KalibrRadtanReprojectionError(double observed_x, double observed_y)
		: observed_x(observed_x), observed_y(observed_y) {
	}

	template <typename T>
	bool operator()(const T* const camera,
		const T* const point,
		T* residuals) const {
		const T* eax = camera;
		const T* C = camera + 3;
		const T* K = camera + 6;
		const T* k = camera + 12;
		const T* r = camera + 17;

		T XC[3];

		XC[0] = point[0] - C[0];
		XC[1] = point[1] - C[1];
		XC[2] = point[2] - C[2];

		T YY[3];
		ceres::AngleAxisRotatePoint(eax, XC, YY);
		//some weird CS change
		T Y[3];
		Y[0] = K[0] * YY[0] + K[1] * YY[1] + K[2] * YY[2];
		Y[1] = K[3] * YY[1] + K[4] * YY[2];
		Y[2] = K[5] * YY[2];

		T rz = T(1) / (Y[2] + k[0] * ceres::sqrt(Y[0] * Y[0] + Y[1] * Y[1] + Y[2] * Y[2]));
		T xz = Y[0] * rz;
		T yz = Y[1] * rz;

		//distortion
		T n = xz*xz + yz*yz;
		T d = n*(r[0] + r[1]*n);
		T p[2];
		p[0] = xz + xz*d + T(2) * r[2] * xz*yz + r[3] * (T(3) * xz*xz + yz*yz);
		p[1] = yz + yz*d + T(2) * r[3] * xz*yz + r[2] * (xz*xz + T(3) * yz*yz);

		T u[2];

		u[0] = k[1] * p[0] + k[3];
		u[1] = k[2] * p[1] + k[4];

		// The error is the difference between the predicted and observed position.
		residuals[0] = u[0] - T(observed_x);
		residuals[1] = u[1] - T(observed_y);

		return true;
	}

	double observed_x;
	double observed_y;
};

#endif