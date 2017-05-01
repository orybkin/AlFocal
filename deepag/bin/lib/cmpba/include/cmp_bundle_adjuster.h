// CMP Bundle Adjustment based on Google Ceres

#pragma once

#include <string>

#ifdef CMP_BUNDLE_ADJUSTMENT_DLLIMPORT
#define DLLAPI __declspec(dllimport)
#else
#define DLLAPI __declspec(dllexport)
#endif


class DLLAPI BAProblemInterface {
public:
	BAProblemInterface();
	double * parameters;
	double * observations;
	double * camera_constraints;
	double * camera_constraint_weights;
	double * error_weights;
	int* point_index_;
	int* camera_index_;
	int n_ref_cams;
	int n_ref_pts;
	int* ref_cam_index;
	double* ref_cam_parameters;
	double* ref_cam_weights;
	int* ref_pts_index;
	double* ref_pts_parameters;
	double* ref_pts_weights;
	int ncams;
	int npts;
	int nobs;
	int proj_func;
	int * fixed_parmask;
	int * shared_parmask;
	int * covariance_mask;
	int num_iterations;
	int num_threads;
	double f_tolerance;
	double g_tolerance;
	double p_tolerance;
	bool robustify;
	bool * constant_cameras;
	bool * constant_points;
	bool constrained_cameras;
	bool free_data;
	bool verbose;
	bool rolling_shutter;
	int * rolling_shutter_orientation; // 0 for vertical, 1 for horizontal
	int * im_size; // image dimensions for calculating the rolling shutter interpolation
	double rolling_shutter_frame_time;
	double rolling_shutter_capture_time;
	double * rolling_shutter_v;
	double * rolling_shutter_w;
	double * rolling_shutter_w_constraints;
	double * rolling_shutter_v_constraints;
	double * rolling_shutter_v_weight;
	double * rolling_shutter_w_weight;
	//quaterionns representing Rv formed from upvector
	double * qRv;
	//upvector Y rotation component
	double *q;

};




DLLAPI typedef struct{	
		std::string input_cams; //Input camera parameter File name//;
		std::string input_pts; ////, //Input points parameters File name//;
		int proj_func; //type of projection function, Options are: 0-snavely, 1-lourakis//;
		std::string parmask; ////, //shared parameter description file name//;
		std::string ref_points; ////, //reference 3D points File name//;
		std::string ref_cams; ////, //reference 3D points File name//;
		std::string camnames; ////, //names of images in alphabetical order, used for printing P matrices//;
		bool estimate_radial;//estimation of radial distortion//;
		bool print_sba; //print output in sba format - THE SAME FORMAT AS INPUT
		bool print_matrices;// print P matrices of cameras - EXPERIMENTAL//;
		bool print_bundler;// print bundle.out//;
		bool print_wrl;// print a wrl file - EXPERIMENTAL//;
		std::string output_cams; //output cameras file name
		std::string output_pts; //output points file name
		std::string output_bundler; //bundler v0.3 format output file name
		std::string output_wrl; //wrl format output file name
		bool constant_cameras; //keep camera parameters constant, adjust only 3D points
		bool constant_points; //keep the 3D structure constant and adjust only cameras

		//ceres original
		std::string trust_region_strategy; //levenberg_marquardt//,
					  //Options are: levenberg_marquardt, dogleg.//;
		std::string dogleg; //traditional_dogleg//, //Options are: traditional_dogleg,//
					  //subspace_dogleg.//;

		bool inner_iterations; //Use inner iterations to non-linearly //
					//refine each successful trust region step.//;

		std::string blocks_for_inner_iterations; //automatic//, //Options are: //
					//automatic, cameras, points, cameras,points, points,cameras//;

		std::string linear_solver; //sparse_schur//, //Options are: //
					  //sparse_schur, dense_schur, iterative_schur, sparse_normal_cholesky, //
					  //dense_qr, dense_normal_cholesky and cgnr.//;
		std::string preconditioner; //jacobi//, //Options are: //
					  //identity, jacobi, schur_jacobi, cluster_jacobi, //
					  //cluster_tridiagonal.//;
		std::string sparse_linear_algebra_library; //suite_sparse//,
					  //Options are: suite_sparse and cx_sparse.//;
		std::string ordering; //automatic//, //Options are: automatic, user.//;

		bool robustify; //Use a robust loss function.//;

		double eta; //Default value for eta. Eta determines the //
					 //accuracy of each linear solve of the truncated newton step. //
					 //Changing this parameter can affect solve performance.//;

		bool use_block_amd; //Use a block oriented fill reducing //
					//ordering.//;

		int num_threads; //Number of threads.//;
		int num_iterations; //Number of iterations.//;
		double max_solver_time; //Maximum solve time in seconds.//;
		bool nonmonotonic_steps; //Trust region algorithm can use//
		
		std::string solver_log; ////, //File to record the solver execution to.//;

}BundleAdjustmentSettings;



DLLAPI void initBundleAdjustmentSettings(BundleAdjustmentSettings * settings);
DLLAPI void processBundleAdjustmentSettings(BundleAdjustmentSettings * settings);
DLLAPI int run_bundle_adjustment();
DLLAPI int SolveProblemFromInterface(BAProblemInterface * problemInterface, std::string * summary = NULL);
