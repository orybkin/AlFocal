#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>



#include "cmp_bundle_adjuster.h"
#include "bundle_adjuster.h"
#include "ceres/ceres.h"
#include "iostream"
#include "fstream"

#include <mex.h>

namespace ceres {
	namespace examples {






		BAProblemInterface * processMatlabInterface(const mxArray * matlab_ceres_interface){
			int n;
			double *pr;
			mxArray * mx;
			BAProblemInterface * problemInterface = new BAProblemInterface();

			mx = mxGetField(matlab_ceres_interface, 0, "verbose");
			if (mx == NULL) {
				problemInterface->verbose = false;
			}
			else{
				problemInterface->verbose = (bool)(*mxGetPr(mx));
				if (problemInterface->verbose)mexPrintf("verbose = %d\n", problemInterface->verbose);
			}

			mx = mxGetField(matlab_ceres_interface, 0, "parameters");
			if (mx == NULL) {
				mexErrMsgTxt("processMatlabInterface(): interface.parameters not found");
				return NULL;
			}
			n = mxGetNumberOfElements(mx);
			problemInterface->parameters = new double[n];
			pr = mxGetPr(mx);
			for (int i = 0; i < n; i++){
				problemInterface->parameters[i] = pr[i];
			}

			mx = mxGetField(matlab_ceres_interface, 0, "camera_constraints");
			if (mx == NULL) {
				problemInterface->shared_parmask = NULL;
			}
			else{
				n = mxGetNumberOfElements(mx);
				problemInterface->camera_constraints = new double[n];
				pr = mxGetPr(mx);
				for (int i = 0; i < n; i++){
					problemInterface->camera_constraints[i] = pr[i];
				}
			}

			mx = mxGetField(matlab_ceres_interface, 0, "camera_constraint_weights");
			if (mx == NULL) {
				problemInterface->shared_parmask = NULL;
			}
			else{
				n = mxGetNumberOfElements(mx);
				problemInterface->camera_constraint_weights = new double[n];
				pr = mxGetPr(mx);
				for (int i = 0; i < n; i++){
					problemInterface->camera_constraint_weights[i] = pr[i];
				}
			}

			mx = mxGetField(matlab_ceres_interface, 0, "observations");
			if (mx == NULL) {
				mexErrMsgTxt("processMatlabInterface(): interface.observations not found");
				return NULL;
			}
			n = mxGetNumberOfElements(mx);
			problemInterface->observations = new double[n];
			pr = mxGetPr(mx);
			for (int i = 0; i < n; i++){
				problemInterface->observations[i] = pr[i];
			}

			mx = mxGetField(matlab_ceres_interface, 0, "camera_index");
			if (mx == NULL) {
				mexErrMsgTxt("processMatlabInterface(): interface.camera_index not found");
				return NULL;
			}
			n = mxGetNumberOfElements(mx);
			problemInterface->camera_index_ = new int[n];
			pr = mxGetPr(mx);
			for (int i = 0; i < n; i++){
				problemInterface->camera_index_[i] = (int)pr[i];
			}

			mx = mxGetField(matlab_ceres_interface, 0, "point_index");
			if (mx == NULL) {
				mexErrMsgTxt("processMatlabInterface(): interface.point_index not found");
				return NULL;
			}
			n = mxGetNumberOfElements(mx);
			problemInterface->point_index_ = new int[n];
			pr = mxGetPr(mx);
			for (int i = 0; i < n; i++){
				problemInterface->point_index_[i] = (int)pr[i];
			}

			mx = mxGetField(matlab_ceres_interface, 0, "fixed_parmask");
			if (mx == NULL) {
				problemInterface->fixed_parmask = NULL;
			}
			else{
				n = mxGetNumberOfElements(mx);
				problemInterface->fixed_parmask = new int[n];
				pr = mxGetPr(mx);
				for (int i = 0; i < n; i++){
					problemInterface->fixed_parmask[i] = (int)pr[i];
				}
			}

			mx = mxGetField(matlab_ceres_interface, 0, "shared_parmask");
			if (mx == NULL) {
				problemInterface->shared_parmask = NULL;
			}
			else{
				n = mxGetNumberOfElements(mx);
				problemInterface->shared_parmask = new int[n];
				pr = mxGetPr(mx);
				for (int i = 0; i < n; i++){
					problemInterface->shared_parmask[i] = (int)pr[i];
				}
			}

			mx = mxGetField(matlab_ceres_interface, 0, "covariance_mask");
			if (mx == NULL) {
				problemInterface->covariance_mask = NULL;
			}
			else{
				n = mxGetNumberOfElements(mx);
				problemInterface->covariance_mask = new int[n];
				pr = mxGetPr(mx);
				for (int i = 0; i < n; i++){
					problemInterface->covariance_mask[i] = (int)pr[i];
				}
			}

			mx = mxGetField(matlab_ceres_interface, 0, "proj_func");
			if (mx == NULL) {
				if (problemInterface->verbose)mexPrintf("proj_func not found, using proj_func=1\n");
				problemInterface->proj_func = 1;

			}
			else{
				problemInterface->proj_func = (int)(*mxGetPr(mx));
			}

			mx = mxGetField(matlab_ceres_interface, 0, "nobs");
			if (mx == NULL) {
				mexErrMsgTxt("processMatlabInterface(): interface.nobs not found\n");
				return NULL;
			}
			problemInterface->nobs = (int)(*mxGetPr(mx));

			mx = mxGetField(matlab_ceres_interface, 0, "ncams");
			if (mx == NULL) {
				mexErrMsgTxt("processMatlabInterface(): interface.ncams not found\n");
				return NULL;
			}
			problemInterface->ncams = (int)(*mxGetPr(mx));

			mx = mxGetField(matlab_ceres_interface, 0, "npts");
			if (mx == NULL) {
				mexErrMsgTxt("processMatlabInterface(): interface.npts not found\n");
				return NULL;
			}
			problemInterface->npts = (int)(*mxGetPr(mx));

			mx = mxGetField(matlab_ceres_interface, 0, "num_iterations");
			if (mx == NULL) {
				problemInterface->num_iterations = 5;
				if (problemInterface->verbose)mexPrintf("max_num_iterations not specified - using 5\n");
			}
			else{
				problemInterface->num_iterations = (int)(*mxGetPr(mx));
			}

			mx = mxGetField(matlab_ceres_interface, 0, "num_threads");
			if (mx == NULL) {
				problemInterface->num_threads = 1;
			}
			else{
				problemInterface->num_threads = (int)(*mxGetPr(mx));
			}
			if (problemInterface->verbose)mexPrintf("%d threads\n", problemInterface->num_threads);

			mx = mxGetField(matlab_ceres_interface, 0, "robustify");
			if (mx == NULL) {
				problemInterface->robustify = 0;
			}
			else{
				problemInterface->robustify = (int)(*mxGetPr(mx));
				if (problemInterface->verbose)mexPrintf("robustify -> %d\n", problemInterface->robustify);
			}
			mx = mxGetField(matlab_ceres_interface, 0, "p_tolerance");
			if (mx == NULL) {
				problemInterface->p_tolerance = 0.001;
			}
			else{
				problemInterface->p_tolerance = (double)(*mxGetPr(mx));
				if (problemInterface->verbose)mexPrintf("p_tolerance = %f\n", problemInterface->p_tolerance);
			}
			mx = mxGetField(matlab_ceres_interface, 0, "g_tolerance");
			if (mx == NULL) {
				problemInterface->g_tolerance = 0.001;
			}
			else{
				problemInterface->g_tolerance = (double)(*mxGetPr(mx));
				if (problemInterface->verbose)mexPrintf("g_tolerance = %f\n", problemInterface->g_tolerance);
			}
			mx = mxGetField(matlab_ceres_interface, 0, "f_tolerance");
			if (mx == NULL) {
				problemInterface->f_tolerance = 0.001;
			}
			else{
				problemInterface->f_tolerance = (double)(*mxGetPr(mx));
				if (problemInterface->verbose)mexPrintf("f_tolerance = %f\n", problemInterface->f_tolerance);
			}
			mx = mxGetField(matlab_ceres_interface, 0, "rolling_shutter");
			if (mx == NULL) {
				problemInterface->rolling_shutter = false;
			}
			else{
				problemInterface->rolling_shutter = (int)(*mxGetPr(mx));
				if (problemInterface->verbose)mexPrintf("rolling_shutter = %d\n", problemInterface->rolling_shutter);
			}
			mx = mxGetField(matlab_ceres_interface, 0, "constant_cameras");
			if (mx == NULL) {
				problemInterface->constant_cameras = NULL;
			}
			else{

				problemInterface->constant_cameras = new bool[problemInterface->ncams];
				pr = mxGetPr(mx);
				for (int i = 0; i < problemInterface->ncams; i++){
					problemInterface->constant_cameras[i] = (bool)pr[i];
				}
				//mexPrintf("constant_cameras = %d\n", problemInterface->constant_cameras);
			}
			mx = mxGetField(matlab_ceres_interface, 0, "constant_points");
			if (mx == NULL) {
				problemInterface->constant_points = NULL;
			}
			else{
				problemInterface->constant_points = new bool[problemInterface->npts];
				pr = mxGetPr(mx);
				for (int i = 0; i < problemInterface->npts; i++){
					problemInterface->constant_points[i] = (bool)pr[i];
				}
			}
			mx = mxGetField(matlab_ceres_interface, 0, "constrained_cameras");
			if (mx == NULL) {
				problemInterface->constrained_cameras = false;
			}
			else{
				problemInterface->constrained_cameras = (int)(*mxGetPr(mx));
				if (problemInterface->verbose)mexPrintf("constrained_cameras = %f\n", problemInterface->constrained_cameras);
			}
			mx = mxGetField(matlab_ceres_interface, 0, "free_data");
			if (mx == NULL) {
				problemInterface->free_data = false;
			}
			else{
				problemInterface->free_data = (int)(*mxGetPr(mx));
				if (problemInterface->verbose)mexPrintf("free_data = %f\n", problemInterface->free_data);
			}

			mx = mxGetField(matlab_ceres_interface, 0, "rolling_shutter_orientation");
			if (mx == NULL) {
				problemInterface->rolling_shutter_orientation = NULL;
			}
			else{
				n = mxGetNumberOfElements(mx);
				problemInterface->rolling_shutter_orientation = new int[n];
				pr = mxGetPr(mx);
				for (int i = 0; i < n; i++){
					problemInterface->rolling_shutter_orientation[i] = (int)pr[i];
				}
			}

			mx = mxGetField(matlab_ceres_interface, 0, "im_size");
			if (mx == NULL) {
				problemInterface->im_size = NULL;
			}
			else{
				n = mxGetNumberOfElements(mx);
				problemInterface->im_size = new int[n];
				pr = mxGetPr(mx);
				for (int i = 0; i < n; i++){
					problemInterface->im_size[i] = (int)pr[i];
				}
			}
			mx = mxGetField(matlab_ceres_interface, 0, "rolling_shutter_frame_time");
			if (mx == NULL) {
				problemInterface->rolling_shutter_frame_time = false;
			}
			else{
				problemInterface->rolling_shutter_frame_time = (double)(*mxGetPr(mx));
				if (problemInterface->verbose)mexPrintf("rolling_shutter_frame_time = %f\n", problemInterface->rolling_shutter_frame_time);
			}

			mx = mxGetField(matlab_ceres_interface, 0, "rolling_shutter_capture_time");
			if (mx == NULL) {
				problemInterface->rolling_shutter_capture_time = false;
			}
			else{
				problemInterface->rolling_shutter_capture_time = (double)(*mxGetPr(mx));
				if (problemInterface->verbose)mexPrintf("rolling_shutter_capture_time = %f\n", problemInterface->rolling_shutter_capture_time);
			}

			mx = mxGetField(matlab_ceres_interface, 0, "rolling_shutter_v");
			if (mx == NULL) {
				problemInterface->rolling_shutter_v = NULL;
			}
			else{

				problemInterface->rolling_shutter_v = new double[problemInterface->ncams * 3];
				pr = mxGetPr(mx);
				for (int i = 0; i < problemInterface->ncams * 3; i++){
					problemInterface->rolling_shutter_v[i] = (double)pr[i];
				}
			}
			mx = mxGetField(matlab_ceres_interface, 0, "rolling_shutter_w");
			if (mx == NULL) {
				problemInterface->rolling_shutter_w = NULL;
			}
			else{

				problemInterface->rolling_shutter_w = new double[problemInterface->ncams * 3];
				pr = mxGetPr(mx);
				for (int i = 0; i < problemInterface->ncams * 3; i++){
					problemInterface->rolling_shutter_w[i] = (double)pr[i];
				}
			}
			mx = mxGetField(matlab_ceres_interface, 0, "rolling_shutter_v_constraints");
			if (mx == NULL) {
				problemInterface->rolling_shutter_v_constraints = NULL;
			}
			else{

				problemInterface->rolling_shutter_v_constraints = new double[problemInterface->ncams * 3];
				pr = mxGetPr(mx);
				for (int i = 0; i < problemInterface->ncams * 3; i++){
					problemInterface->rolling_shutter_v_constraints[i] = (double)pr[i];
				}
			}

			mx = mxGetField(matlab_ceres_interface, 0, "rolling_shutter_w_constraints");
			if (mx == NULL) {
				problemInterface->rolling_shutter_w_constraints = NULL;
			}
			else{

				problemInterface->rolling_shutter_w_constraints = new double[problemInterface->ncams * 3];
				pr = mxGetPr(mx);
				for (int i = 0; i < problemInterface->ncams * 3; i++){
					problemInterface->rolling_shutter_w_constraints[i] = (double)pr[i];
				}
			}

			mx = mxGetField(matlab_ceres_interface, 0, "rolling_shutter_v_weight");
			if (mx == NULL) {
				problemInterface->rolling_shutter_v_weight = NULL;
			}
			else{

				problemInterface->rolling_shutter_v_weight = new double[problemInterface->ncams * 3];
				pr = mxGetPr(mx);
				for (int i = 0; i < problemInterface->ncams * 3; i++){
					problemInterface->rolling_shutter_v_weight[i] = (double)pr[i];
				}
			}

			mx = mxGetField(matlab_ceres_interface, 0, "rolling_shutter_w_weight");
			if (mx == NULL) {
				problemInterface->rolling_shutter_w_weight = NULL;
			}
			else{

				problemInterface->rolling_shutter_w_weight = new double[problemInterface->ncams * 3];
				pr = mxGetPr(mx);
				for (int i = 0; i < problemInterface->ncams * 3; i++){
					problemInterface->rolling_shutter_w_weight[i] = (double)pr[i];
				}
			}

			mx = mxGetField(matlab_ceres_interface, 0, "n_ref_cams");
			if (mx == NULL) {
				problemInterface->n_ref_cams = 0;
			}
			else{
				problemInterface->n_ref_cams = (int)(*mxGetPr(mx));
			}	

			mx = mxGetField(matlab_ceres_interface, 0, "n_ref_pts");
			if (mx == NULL) {
				problemInterface->n_ref_pts = 0;
			}
			else{
				problemInterface->n_ref_pts = (int)(*mxGetPr(mx));
			}

			mx = mxGetField(matlab_ceres_interface, 0, "ref_cam_parameters");
			if (mx == NULL) {
				problemInterface->ref_cam_parameters = NULL;
			}
			else{
				problemInterface->ref_cam_parameters = new double[problemInterface->n_ref_cams*3];
				pr = mxGetPr(mx);
				for (int i = 0; i < problemInterface->n_ref_cams*3; i++){
					problemInterface->ref_cam_parameters[i] = (double)pr[i];
				}
			}

			mx = mxGetField(matlab_ceres_interface, 0, "ref_pts_parameters");
			if (mx == NULL) {
				problemInterface->ref_pts_parameters = NULL;
			}
			else{
				problemInterface->ref_pts_parameters = new double[problemInterface->n_ref_pts * 3];
				pr = mxGetPr(mx);
				for (int i = 0; i < problemInterface->n_ref_pts * 3; i++){
					problemInterface->ref_pts_parameters[i] = (double)pr[i];
				}
			}

			mx = mxGetField(matlab_ceres_interface, 0, "ref_cam_index");
			if (mx == NULL) {
				problemInterface->ref_cam_index = NULL;
			}
			else{
				problemInterface->ref_cam_index = new int[problemInterface->n_ref_cams];
				pr = mxGetPr(mx);
				for (int i = 0; i < problemInterface->n_ref_cams; i++){
					problemInterface->ref_cam_index[i] = (int)pr[i];
				}
			}
			
			mx = mxGetField(matlab_ceres_interface, 0, "ref_pts_index");
			if (mx == NULL) {
				problemInterface->ref_pts_index = NULL;
			}
			else{
				problemInterface->ref_pts_index = new int[problemInterface->n_ref_pts];
				pr = mxGetPr(mx);
				for (int i = 0; i < problemInterface->n_ref_pts; i++){
					problemInterface->ref_pts_index[i] = (int)pr[i];
				}
			}

			mx = mxGetField(matlab_ceres_interface, 0, "ref_cam_weights");
			if (mx == NULL) {
				problemInterface->ref_cam_weights = NULL;
			}
			else{
				problemInterface->ref_cam_weights = new double[problemInterface->n_ref_cams];
				pr = mxGetPr(mx);
				for (int i = 0; i < problemInterface->n_ref_cams; i++){
					problemInterface->ref_cam_weights[i] = (double)pr[i];
				}
			}

			mx = mxGetField(matlab_ceres_interface, 0, "ref_pts_weights");
			if (mx == NULL) {
				problemInterface->ref_pts_weights = NULL;
			}
			else{
				problemInterface->ref_pts_weights = new double[problemInterface->n_ref_pts];
				pr = mxGetPr(mx);
				for (int i = 0; i < problemInterface->n_ref_pts; i++){
					problemInterface->ref_pts_weights[i] = (double)pr[i];
				}
			}

			mx = mxGetField(matlab_ceres_interface, 0, "error_weigths");
			if (mx == NULL) {
				problemInterface->error_weights = NULL;
			}
			else{
				problemInterface->error_weights = new double[problemInterface->nobs];
				pr = mxGetPr(mx);
				for (int i = 0; i < problemInterface->nobs; i++){
					problemInterface->error_weights[i] = (double)pr[i];
				}
			}



			return problemInterface;
		}


		int SolveProblemMex(const mxArray * matlab_ceres_interface, int nlhs, mxArray * plhs[]) {

			if (!mxIsStruct(matlab_ceres_interface)) {
				mexErrMsgTxt("SolveProblem(): interface parameter must be struct");
				return 1;
			}



			BAProblem ba_problem;
			BAProblemInterface * problemInterface;
			problemInterface = processMatlabInterface(matlab_ceres_interface);
			if (problemInterface == NULL){
				mexErrMsgTxt("SolveProblem(): failed to process the interface");
				return 1;
			}
			Problem problem;
			ba_problem.loadParametersFromInterface(problemInterface);



			Solver::Options options;
			options.function_tolerance = problemInterface->f_tolerance;
			options.parameter_tolerance = problemInterface->p_tolerance;
			options.gradient_tolerance = problemInterface->g_tolerance;
			options.max_num_iterations = problemInterface->num_iterations;
			options.num_threads = problemInterface->num_threads;
			Solver::Summary summary;
			std::string report;
			SolveProblem(&ba_problem, options,&summary);
			if (problemInterface->verbose){
				mexPrintf(summary.FullReport().c_str());
				mexPrintf("\n");
			}

			//std::cout << summary.FullReport() << "\n";
			//return parameters
			if (nlhs > 0){
				plhs[0] = mxCreateNumericMatrix(1, ba_problem.num_parameters_, mxDOUBLE_CLASS, mxREAL);
				double * data = (double *)mxGetData(plhs[0]);
				for (int i = 0; i < ba_problem.num_parameters_; i++){
					data[i] = ba_problem.parameters_[i];
				}
			}
			if (nlhs > 1){
				plhs[1] = mxCreateDoubleMatrix(summary.iterations.size(), 1, mxREAL);
				double * data = (double *)mxGetData(plhs[1]);
				for (int i = 0; i < summary.iterations.size(); i++){
					data[i] = summary.iterations[i].cost;
				}

			}
			

			//RS parameters
			if (ba_problem.rolling_shutter == 1){
				if (nlhs > 2){
					plhs[2] = mxCreateNumericMatrix(1, ba_problem.num_cameras_*3, mxDOUBLE_CLASS, mxREAL);
					double * data = (double *)mxGetData(plhs[2]);
					for (int i = 0; i < ba_problem.num_cameras_ * 3; i++){
						data[i] = ba_problem.rolling_shutter_v[i];
					}
				}
				if (nlhs > 3){
					plhs[3] = mxCreateDoubleMatrix(1, ba_problem.num_cameras_*3, mxREAL);
					double * data = (double *)mxGetData(plhs[3]);
					for (int i = 0; i < ba_problem.num_cameras_ * 3; i++){
						data[i] = ba_problem.rolling_shutter_w[i];
					}

				}
			
			
			}




			delete problemInterface;
			return 0;



		}

		

	}
};

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]){
	if (nrhs < 1){
		mexErrMsgTxt("cmp_bundle_adjuster_mex: not enough input arguments");
	}
	else{
		ceres::examples::SolveProblemMex(prhs[0], nlhs, plhs);
	}



}