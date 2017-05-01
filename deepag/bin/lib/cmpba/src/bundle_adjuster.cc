// Ceres Solver - A fast non-linear least squares minimizer
// Copyright 2010, 2011, 2012 Google Inc. All rights reserved.
// http://code.google.com/p/ceres-solver/
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the following disclaimer.
// * Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
// * Neither the name of Google Inc. nor the names of its contributors may be
//   used to endorse or promote products derived from this software without
//   specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// Author: sameeragarwal@google.com (Sameer Agarwal)
//
// An example of solving a dynamically sized problem with various
// solvers and loss functions.
//
// For a simpler bare bones example of doing bundle adjustment with
// Ceres, please see simple_bundle_adjuster.cc.
//
// NOTE: This example will not compile without gflags and SuiteSparse.
//
// The problem being solved here is known as a Bundle Adjustment
// problem in computer vision. Given a set of 3d points X_1, ..., X_n,
// a set of cameras P_1, ..., P_m. If the point X_i is visible in
// image j, then there is a 2D observation u_ij that is the expected
// projection of X_i using P_j. The aim of this optimization is to
// find values of X_i and P_j such that the reprojection error
//
//    E(X,P) =  sum_ij  |u_ij - P_j X_i|^2
//
// is minimized.
//
// The problem used here comes from a collection of bundle adjustment
// problems published at University of Washington.
// http://grail.cs.washington.edu/projects/bal

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "bundle_adjuster.h"


using namespace ceres;
using namespace std;

//cmp added
DEFINE_string(input_cams, "", "Input camera parameter File name");
DEFINE_string(input_pts, "", "Input points parameters File name");
DEFINE_string(input_boards, "", "Input parameters of calibration boards File name");
DEFINE_string(input_boards_pts, "", "Input parameters of calibration boards File name");
DEFINE_int32(proj_func, 1, "type of projection function, Options are: 0-snavely, 1-lourakis");
DEFINE_string(parmask, "", "shared parameter description file name");
DEFINE_string(ref_points, "", "reference 3D points File name");
DEFINE_string(ref_cams, "", "reference 3D cameras File name");
DEFINE_string(camnames, "", "names of images in alphabetical order, used for printing P matrices");
DEFINE_bool(estimate_radial,false,"estimation of radial distortion");
DEFINE_bool(print_matrices,false," print P matrices of cameras");
DEFINE_bool(print_bundler,false," print output in bundler v0.3 format");
DEFINE_bool(print_wrl,false," print output in wrl format");
DEFINE_bool(print_sba,false," print output in sba format");
DEFINE_bool(print_calib,false," print output in calib format");
DEFINE_string(output_cams,"cams_res.txt", "output cameras file name");
DEFINE_string(output_pts,"pts_res.txt", "output points file name");
DEFINE_string(output_boards,"boards_res.txt", "calibration boards parameters");
DEFINE_string(output_bundler,"bundle.out", "bundler v0.3 format output file name");
DEFINE_string(output_wrl,"bundle.wrl", "wrl format output file name");
DEFINE_bool(constant_cameras,false,"keep the camera positions constant");
DEFINE_bool(constant_points,false,"keep the 3D structure constant and adjust only cameras");
DEFINE_int32(fixed_par, 0, "fix parameters?");
DEFINE_int32(rolling_shutter, 0, "rolling shutter?");
DEFINE_int32(rs_direction, 0, "0-vertical,1-horizontal");



DEFINE_int32(dump,0, "0 - don't print iterations, 1 - print them");

//ceres original
DEFINE_string(input, "", "Input File name");
DEFINE_string(trust_region_strategy, "levenberg_marquardt",
	"Options are: levenberg_marquardt, dogleg.");
DEFINE_string(dogleg, "traditional_dogleg", "Options are: traditional_dogleg,"
	"subspace_dogleg.");

DEFINE_bool(inner_iterations, false, "Use inner iterations to non-linearly "
	"refine each successful trust region step.");

DEFINE_string(blocks_for_inner_iterations, "automatic", "Options are: "
	"automatic, cameras, points, cameras,points, points,cameras");

DEFINE_string(linear_solver, "sparse_schur", "Options are: "
	"sparse_schur, dense_schur, iterative_schur, sparse_normal_cholesky, "
	"dense_qr, dense_normal_cholesky and cgnr.");
DEFINE_string(preconditioner, "jacobi", "Options are: "
	"identity, jacobi, schur_jacobi, cluster_jacobi, "
	"cluster_tridiagonal.");
DEFINE_string(sparse_linear_algebra_library, "suite_sparse",
	"Options are: suite_sparse and cx_sparse.");
DEFINE_string(ordering, "automatic", "Options are: automatic, user.");

DEFINE_bool(use_quaternions, false, "If true, uses quaternions to represent "
	"rotations. If false, angle axis is used.");
DEFINE_bool(use_local_parameterization, false, "For quaternions, use a local "
	"parameterization.");
DEFINE_bool(robustify, false, "Use a robust loss function.");

DEFINE_double(eta, 1e-2, "Default value for eta. Eta determines the "
	"accuracy of each linear solve of the truncated newton step. "
	"Changing this parameter can affect solve performance.");

DEFINE_bool(use_block_amd, true, "Use a block oriented fill reducing "
	"ordering.");

DEFINE_int32(num_threads, 1, "Number of threads.");
DEFINE_int32(num_iterations, 5, "Number of iterations.");
DEFINE_double(max_solver_time, 1e32, "Maximum solve time in seconds.");
DEFINE_bool(nonmonotonic_steps, false, "Trust region algorithm can use"
	" nonmonotic steps.");

DEFINE_double(rotation_sigma, 0.0, "Standard deviation of camera rotation "
	"perturbation.");
DEFINE_double(translation_sigma, 0.0, "Standard deviation of the camera "
	"translation perturbation.");
DEFINE_double(point_sigma, 0.0, "Standard deviation of the point "
	"perturbation.");
DEFINE_int32(random_seed, 38401, "Random seed used to set the state "
	"of the pseudo random number generator used to generate "
	"the pertubations.");
DEFINE_string(solver_log, "log.txt", "File to record the solver execution to.");




		void SetLinearSolver(Solver::Options* options) {
			CHECK(StringToLinearSolverType(FLAGS_linear_solver,
				&options->linear_solver_type));
			CHECK(StringToPreconditionerType(FLAGS_preconditioner,
				&options->preconditioner_type));
			CHECK(StringToSparseLinearAlgebraLibraryType(
				FLAGS_sparse_linear_algebra_library,
				&options->sparse_linear_algebra_library_type));
			options->num_linear_solver_threads = FLAGS_num_threads;
			
		}


		void SetOrdering(BAProblem* bal_problem, Solver::Options* options) {
			const int num_points = bal_problem->num_points_;
			const int point_block_size = bal_problem->num_point_parameters_;
			double* points = bal_problem->mutable_points();

			const int num_cameras = bal_problem->num_cameras_;
			const int camera_block_size = bal_problem->num_point_parameters_;
			double* cameras = bal_problem->mutable_cameras();


			if (options->use_inner_iterations) {
				if (FLAGS_blocks_for_inner_iterations == "cameras") {
					LOG(INFO) << "Camera blocks for inner iterations";
					options->inner_iteration_ordering.reset(new ParameterBlockOrdering);
					for (int i = 0; i < num_cameras; ++i) {
						options->inner_iteration_ordering->AddElementToGroup(cameras + camera_block_size * i, 0);
					}
				} else if (FLAGS_blocks_for_inner_iterations == "points") {
					LOG(INFO) << "Point blocks for inner iterations";
					options->inner_iteration_ordering.reset(new ParameterBlockOrdering);
					for (int i = 0; i < num_points; ++i) {
						options->inner_iteration_ordering->AddElementToGroup(points + point_block_size * i, 0);
					}
				} else if (FLAGS_blocks_for_inner_iterations == "cameras,points") {
					LOG(INFO) << "Camera followed by point blocks for inner iterations";
					options->inner_iteration_ordering.reset(new ParameterBlockOrdering);
					for (int i = 0; i < num_cameras; ++i) {
						options->inner_iteration_ordering->AddElementToGroup(cameras + camera_block_size * i, 0);
					}
					for (int i = 0; i < num_points; ++i) {
						options->inner_iteration_ordering->AddElementToGroup(points + point_block_size * i, 1);
					}
				} else if (FLAGS_blocks_for_inner_iterations == "points,cameras") {
					LOG(INFO) << "Point followed by camera blocks for inner iterations";
					options->inner_iteration_ordering.reset(new ParameterBlockOrdering);
					for (int i = 0; i < num_cameras; ++i) {
						options->inner_iteration_ordering->AddElementToGroup(cameras + camera_block_size * i, 1);
					}
					for (int i = 0; i < num_points; ++i) {
						options->inner_iteration_ordering->AddElementToGroup(points + point_block_size * i, 0);
					}
				} else if (FLAGS_blocks_for_inner_iterations == "automatic") {
					LOG(INFO) << "Choosing automatic blocks for inner iterations";
				} else {
					LOG(FATAL) << "Unknown block type for inner iterations: "
						<< FLAGS_blocks_for_inner_iterations;
				}
			}

			// Bundle adjustment problems have a sparsity structure that makes
			// them amenable to more specialized and much more efficient
			// solution strategies. The SPARSE_SCHUR, DENSE_SCHUR and
			// ITERATIVE_SCHUR solvers make use of this specialized
			// structure.
			//
			// This can either be done by specifying Options::ordering_type =
			// ceres::SCHUR, in which case Ceres will automatically determine
			// the right ParameterBlock ordering, or by manually specifying a
			// suitable ordering vector and defining
			// Options::num_eliminate_blocks.
			if (FLAGS_ordering == "automatic") {
				return;
			}

			ceres::ParameterBlockOrdering* ordering =
				new ceres::ParameterBlockOrdering;

			// The points come before the cameras.
			for (int i = 0; i < num_points; ++i) {
				ordering->AddElementToGroup(points + point_block_size * i, 0);
			}

			for (int i = 0; i < num_cameras; ++i) {
				// When using axis-angle, there is a single parameter block for
				// the entire camera.
				ordering->AddElementToGroup(cameras + camera_block_size * i, 1);
				//ordering->AddElementToGroup(cameras + camera_block_size * i+3, 1);
				// If quaternions are used, there are two blocks, so add the
				// second block to the ordering.
				if (FLAGS_use_quaternions) {
					ordering->AddElementToGroup(cameras + camera_block_size * i + 4, 1);
				}
			}

			options->linear_solver_ordering.reset(ordering);
		}

		void SetMinimizerOptions(Solver::Options* options) {
			options->max_num_iterations = FLAGS_num_iterations;
			options->minimizer_progress_to_stdout = true;
			options->num_threads = FLAGS_num_threads;
			options->eta = FLAGS_eta;
			options->max_solver_time_in_seconds = FLAGS_max_solver_time;
			options->use_nonmonotonic_steps = FLAGS_nonmonotonic_steps;
			CHECK(StringToTrustRegionStrategyType(FLAGS_trust_region_strategy,
				&options->trust_region_strategy_type));
			CHECK(StringToDoglegType(FLAGS_dogleg, &options->dogleg_type));
			options->use_inner_iterations = FLAGS_inner_iterations;
		}

		void SetSolverOptionsFromFlags(BAProblem* bal_problem,
			Solver::Options* options) {
				SetMinimizerOptions(options);
				SetLinearSolver(options);
				SetOrdering(bal_problem, options);
		}

		void addReferencePoints(BAProblem * ba_problem, Problem* problem){

			double * points = ba_problem->mutable_points();
			//ba_problem->loadReferencePoints(FLAGS_ref_points);
			double * ref_pts = ba_problem->ref_pts_;
			int * pt_ids = ba_problem->ref_pts_id_;
			double * ref_pts_weights = ba_problem->ref_pts_weights_;

			for (int i = 0; i<ba_problem->n_ref_pts_; i++){
				CostFunction* cost_function;
				cost_function = new AutoDiffCostFunction<Reference3DPointError, 3, 3>(
					new Reference3DPointError(ref_pts + i * 3, ref_pts_weights[i]));

				problem->AddResidualBlock(cost_function,
					NULL,
					points + pt_ids[i] * 3
					);
			}
		}

		void addReferenceCameras(BAProblem * sba_problem, Problem* problem){

			double * cameras = sba_problem->mutable_cameras();
			sba_problem->loadReferenceCameras(FLAGS_ref_cams);
			double * ref_cams = sba_problem->ref_cams_;
			double * ref_cams_weights = sba_problem->ref_cams_weights_;
			int * cams_ids = sba_problem->ref_cams_id_;
			int size = 7;
			for (int i = 0; i<sba_problem->n_ref_cams_; i++){
				CostFunction* cost_function;

				switch (sba_problem->proj_func_){
				case 0:
					cost_function = new AutoDiffCostFunction<ReferenceCamera6pError, 3, 3>(
						new ReferenceCamera6pError(ref_cams + i * 3, ref_cams_weights[i]));
					break;
				case 1:
					cost_function =
						cost_function = new AutoDiffCostFunction<ReferenceCamera6pError, 3, 11>(
						new ReferenceCamera6pError(ref_cams + i * 3, ref_cams_weights[i]));
					break;
				case 2:
					cost_function =
						cost_function = new AutoDiffCostFunction<ReferenceCamera6pError, 3, 6>(
						new ReferenceCamera6pError(ref_cams + i * 3, ref_cams_weights[i]));
					break;
				case 3:
					cost_function =
						cost_function = new AutoDiffCostFunction<ReferenceCamera6pError, 3, 9>(
						new ReferenceCamera6pError(ref_cams + i * 3, ref_cams_weights[i]));
					break;
				case 4:
					cost_function =
						cost_function = new AutoDiffCostFunction<ReferenceCamera6pError, 3, 10>(
						new ReferenceCamera6pError(ref_cams + i * 3, ref_cams_weights[i]));
					break;
				case 5:
					cost_function =
						cost_function = new AutoDiffCostFunction<ReferenceCamera6pError, 3, 7>(
						new ReferenceCamera6pError(ref_cams + i * 3, ref_cams_weights[i]));
					break;
				case 6:
					cost_function =
						cost_function = new AutoDiffCostFunction<ReferenceCamera6pError, 3, 13>(
						new ReferenceCamera6pError(ref_cams + i * 3, ref_cams_weights[i]));
					break;
				case 7:
					cost_function =
						cost_function = new AutoDiffCostFunction<ReferenceCamera6pError, 3, 6>(
						new ReferenceCamera6pError(ref_cams + i * 3, ref_cams_weights[i]));
					break;
				case 8:
					cost_function =
						cost_function = new AutoDiffCostFunction<ReferenceCamera6pError, 3, 11>(
						new ReferenceCamera6pError(ref_cams + i * 3, ref_cams_weights[i]));
					break;
				case 9:
					cost_function =
						cost_function = new AutoDiffCostFunction<ReferenceCamera6pError, 3, 12>(
						new ReferenceCamera6pError(ref_cams + i * 3, ref_cams_weights[i]));
					break;
				case 10:
					cost_function =
						cost_function = new AutoDiffCostFunction<ReferenceCamera6pError, 3, 9>(
						new ReferenceCamera6pError(ref_cams + i * 3, ref_cams_weights[i]));
					break;
				case 20:
					cost_function =
						cost_function = new AutoDiffCostFunction<ReferenceCamera6pError, 3, 7>(
						new ReferenceCamera6pError(ref_cams + i * 3, ref_cams_weights[i]));
					break;
				case 21:
					cost_function =
						cost_function = new AutoDiffCostFunction<ReferenceCamera6pError, 3, 9>(
						new ReferenceCamera6pError(ref_cams + i * 3, ref_cams_weights[i]));
					break;
				case 22:
				case 24:
					cost_function =
						cost_function = new AutoDiffCostFunction<Reference3DPointError, 3, 4>(
						new Reference3DPointError(ref_cams + i * 3, ref_cams_weights[i]));
					break;
				case 23:
				case 25:
					cost_function =
						cost_function = new AutoDiffCostFunction<Reference3DPointError, 3, 6>(
						new Reference3DPointError(ref_cams + i * 3, ref_cams_weights[i]));
					break;
				case 30:
					cost_function =
						cost_function = new AutoDiffCostFunction<ReferenceCamera7pError, 3, 7>(
						new ReferenceCamera7pError(ref_cams + i * 3, ref_cams_weights[i]));
					break;
				case 31:
					cost_function =
						cost_function = new AutoDiffCostFunction<ReferenceCamera7pError, 3, 9>(
						new ReferenceCamera7pError(ref_cams + i * 3, ref_cams_weights[i]));
					break;
				case 40:
					cost_function =
						cost_function = new AutoDiffCostFunction<ReferenceCamera6pError, 3, 9>(
						new ReferenceCamera6pError(ref_cams + i * 3, ref_cams_weights[i]));
					break;
				case 41:
					cost_function =
						cost_function = new AutoDiffCostFunction<ReferenceCamera6pError, 3, 6>(
						new ReferenceCamera6pError(ref_cams + i * 3, ref_cams_weights[i]));
					break;
				case 50:
					cost_function =
						cost_function = new AutoDiffCostFunction<ReferenceCamera6pError, 3, 9>(
						new ReferenceCamera6pError(ref_cams + i * 3, ref_cams_weights[i]));
					break;
				case 51:
					cost_function =
						cost_function = new AutoDiffCostFunction<ReferenceCamera6pError, 3, 6>(
						new ReferenceCamera6pError(ref_cams + i * 3, ref_cams_weights[i]));
					break;
				}


				switch (sba_problem->proj_func_){
				case 0:
					problem->AddResidualBlock(cost_function,
						NULL,
						cameras + cams_ids[i] * sba_problem->num_camera_parameters_ + 4
						);
					break;
				case 1:
				case 2:
				case 3:
				case 4:
				case 5:
				case 6:
				case 7:
				case 8:
				case 9:
				case 10:
				case 20:
				case 21:
					
					problem->AddResidualBlock(cost_function,
						NULL,
						cameras + cams_ids[i] * sba_problem->num_camera_parameters_ + 1
						);
					break;
				case 22:
				case 23:
				case 24:
				case 25:
					problem->AddResidualBlock(cost_function,
						NULL,
						cameras + cams_ids[i] * sba_problem->num_camera_parameters_ + 4
						);
					break;
				case 30:
				case 31:
					problem->AddResidualBlock(cost_function,
						NULL,
						cameras + cams_ids[i] * sba_problem->num_camera_parameters_
						);
					break;
				case 40:
				case 41:
				case 50:
				case 51:
					problem->AddResidualBlock(cost_function,
						NULL,
						cameras + cams_ids[i] * sba_problem->num_camera_parameters_ + 1
						);
					break;
				}



			}



		}


		void BuildProblem(BAProblem* sba_problem, Problem* problem) {
			const int point_block_size = sba_problem->num_point_parameters_;
			const int camera_block_size = sba_problem->num_camera_parameters_;
			double* points = sba_problem->mutable_points();
			double* cameras = sba_problem->mutable_cameras();
			int num_camera_parameters;

			printf("num cameras: %d\n", sba_problem->num_cameras_);
			printf("num points: %d\n", sba_problem->num_points_);
			printf("num observations: %d\n", sba_problem->num_observations_);

			// Observations is 2*num_observations long array observations =
			// [u_1, u_2, ... , u_n], where each u_i is two dimensional, the x
			// and y positions of the observation.
			const double* observations = sba_problem->observations_;
			LocalParameterization* quaternion_parameterization = new QuaternionParameterization;
			
			

			int cnt=0;
			double weight;
			for (int i = 0; i < sba_problem->num_observations_; ++i) {
				CostFunction* cost_function;
				// Each Residual block takes a point and a camera as input and
				// outputs a 2 dimensional residual.
				double* camera = cameras + camera_block_size * sba_problem->camera_index_[i];
				double* next_camera = camera+sba_problem->num_camera_parameters_;
				int nextcam_id = sba_problem->camera_index_[i]+1;

				if(sba_problem->camera_index_[i]==sba_problem->num_cameras_-1){

					next_camera = sba_problem->rs_fake_last_camera;
					nextcam_id = sba_problem->camera_index_[i];
				}

				if (sba_problem->weighted_errors){
					weight = sba_problem->error_weights[i];
				}
				else{
					weight = 1;
				}

				switch(sba_problem->proj_func_){
				case 0:
					cost_function =
						new AutoDiffCostFunction<ReprojectionErrorCustom, 2,3,3,1,2,1,1,2,3>(
						new ReprojectionErrorCustom(observations[2 * i + 0],observations[2 * i + 1],sba_problem->init_rot_+sba_problem->camera_index_[i]*4,weight));
					break;
				case 1:
					num_camera_parameters = 11;
					cost_function =
						new AutoDiffCostFunction<ReprojectionErrorKQT, 2, 11 ,3>(
						new ReprojectionErrorKQT(observations[2 * i + 0], observations[2 * i + 1], sba_problem->init_rot_ + sba_problem->camera_index_[i] * 4, weight));
					break;
				case 2:
					num_camera_parameters = 6;
					cost_function =
						new AutoDiffCostFunction<ReprojectionErrorQT, 2, 6 ,3>(
						new ReprojectionErrorQT(observations[2 * i + 0], observations[2 * i + 1], sba_problem->init_rot_ + sba_problem->camera_index_[i] * 4, camera + 7, weight));
					break;
				case 3:
					num_camera_parameters = 9;
					cost_function =
						new AutoDiffCostFunction<ReprojectionErrorFxX0Y0QT, 2, 9 ,3>(
						new ReprojectionErrorFxX0Y0QT(observations[2 * i + 0], observations[2 * i + 1], sba_problem->init_rot_ + sba_problem->camera_index_[i] * 4, camera[11], camera[10], weight));
					break;
				case 4:
					num_camera_parameters = 10;
					cost_function =
						new AutoDiffCostFunction<ReprojectionErrorArFxX0Y0QT, 2, 10 ,3>(
						new ReprojectionErrorArFxX0Y0QT(observations[2 * i + 0], observations[2 * i + 1], sba_problem->init_rot_ + sba_problem->camera_index_[i] * 4, camera[11], weight));
					break;
				case 5:
					num_camera_parameters = 7;
					cost_function =
						new AutoDiffCostFunction<ReprojectionErrorFQT, 2, 7 ,3>(
						new ReprojectionErrorFQT(observations[2 * i + 0], observations[2 * i + 1], sba_problem->init_rot_ + sba_problem->camera_index_[i] * 4, camera + 8, weight));
					break;
				case 6:
					num_camera_parameters = 13;
					cost_function =
						new AutoDiffCostFunction<ReprojectionErrorKQTRadial, 2,13,3>(
						new ReprojectionErrorKQTRadial(observations[2 * i + 0], observations[2 * i + 1], sba_problem->init_rot_ + sba_problem->camera_index_[i] * 4, weight));
					break;
				case 7:
					num_camera_parameters = 6;
					cost_function =
						new AutoDiffCostFunction<ReprojectionErrorQTRadial, 2,6,3>(
						new ReprojectionErrorQTRadial(observations[2 * i + 0], observations[2 * i + 1], sba_problem->init_rot_ + sba_problem->camera_index_[i] * 4, camera + 9, camera + 7, weight));
					break;
				case 8:
					num_camera_parameters = 11;
					cost_function =
						new AutoDiffCostFunction<ReprojectionErrorFxX0Y0QTRadial, 2,11,3>(
						new ReprojectionErrorFxX0Y0QTRadial(observations[2 * i + 0], observations[2 * i + 1], sba_problem->init_rot_ + sba_problem->camera_index_[i] * 4, camera[13], camera[12], weight));
					break;
				case 9:
					num_camera_parameters = 12;
					cost_function =
						new AutoDiffCostFunction<ReprojectionErrorArFxX0Y0QTRadial, 2,12,3>(
						new ReprojectionErrorArFxX0Y0QTRadial(observations[2 * i + 0], observations[2 * i + 1], sba_problem->init_rot_ + sba_problem->camera_index_[i] * 4, camera[13], weight));
					break;
				case 10:
					num_camera_parameters = 9;
					cost_function =
						new AutoDiffCostFunction<ReprojectionErrorFQTRadial, 2,9,3>(
						new ReprojectionErrorFQTRadial(observations[2 * i + 0], observations[2 * i + 1], sba_problem->init_rot_ + sba_problem->camera_index_[i] * 4, camera + 10, weight));
					break;
				case 11:
					num_camera_parameters = 9;
					cost_function =
						new AutoDiffCostFunction<ReprojectionErrorCalibVisualSfMRadial, 2, 9, 3>(
						new ReprojectionErrorCalibVisualSfMRadial(observations[2 * i + 0], observations[2 * i + 1], sba_problem->init_rot_ + sba_problem->camera_index_[i] * 4));
					break;
				case 12:
					num_camera_parameters = 12;
					cost_function =
						new AutoDiffCostFunction<ReprojectionErrorKannala, 2, 12, 3>(
						new ReprojectionErrorKannala(observations[2 * i + 0], observations[2 * i + 1]));
					break;
				case 20:
					num_camera_parameters = 7;
					cost_function =
						new AutoDiffCostFunction<SnavelyReprojectionError, 2,7,3>(
						new SnavelyReprojectionError(observations[2 * i + 0],observations[2 * i + 1],camera[12],camera[13],camera[10],camera[11]));
					break;
				case 21:
					num_camera_parameters = 9;
					cost_function =
						new AutoDiffCostFunction<SnavelyReprojectionErrorRadial, 2,9,3>(
						new SnavelyReprojectionErrorRadial(observations[2 * i + 0],observations[2 * i + 1],camera[12],camera[13],camera[10],camera[11]));
					break;
				case 22:
				case 24:
					num_camera_parameters = 4;
					cost_function =
						new AutoDiffCostFunction<SnavelyReprojectionErrorWithQuaternions, 2,4,4,3>(
						new SnavelyReprojectionErrorWithQuaternions(observations[2 * i + 0],observations[2 * i + 1],camera[12],camera[13],camera[10],camera[11]));
					break;
				case 23:
				case 25:
					num_camera_parameters = 4;
					cost_function =
						new AutoDiffCostFunction<SnavelyReprojectionErrorWithQuaternionsRadial, 2,4,6,3>(
						new SnavelyReprojectionErrorWithQuaternionsRadial(observations[2 * i + 0],observations[2 * i + 1],camera[12],camera[13],camera[10],camera[11]));
					break;
				case 30: 
					num_camera_parameters = 7;
					cost_function =
						new AutoDiffCostFunction<Quaternion8dReprojectionError, 2,7,3>(
						new Quaternion8dReprojectionError(observations[2 * i + 0],observations[2 * i + 1],camera[10],camera[11],camera[13],camera[9]));
					break;
				case 31: 
					num_camera_parameters = 9;
					cost_function =
						new AutoDiffCostFunction<Quaternion8dReprojectionErrorRadial, 2,9,3>(
						new Quaternion8dReprojectionErrorRadial(observations[2 * i + 0],observations[2 * i + 1],camera[10],camera[11],camera[13],camera[9]));
					break;
				case 40:
					num_camera_parameters = 9;
					// Bundler projection
					cost_function =
						new AutoDiffCostFunction<ReprojectionErrorBundlerRadial, 2,9,3>(
						new ReprojectionErrorBundlerRadial(observations[2 * i + 0],observations[2 * i + 1],sba_problem->init_rot_+sba_problem->camera_index_[i]*4,camera+10));
					break;
				case 41:
				case 42:
					// Bundler projection with shared camera calibration
					num_camera_parameters = 6;
					cost_function =
						new AutoDiffCostFunction<ReprojectionErrorBundlerRadialShared, 2,6,3,3>(
						new ReprojectionErrorBundlerRadialShared(observations[2 * i + 0],observations[2 * i + 1],sba_problem->init_rot_+sba_problem->camera_index_[i]*4,cameras+10));
					break;
				case 50:
					num_camera_parameters = 9;
					if(sba_problem->rolling_shutter_orientation[sba_problem->camera_index_[i]]){
						// Bundler projection with vertical rolling shutter
						cost_function =
							new AutoDiffCostFunction<ReprojectionErrorBundlerRadialRSVideoVertical, 2,9,9,3>(
							new ReprojectionErrorBundlerRadialRSVideoVertical(observations[2 * i + 0],
							observations[2 * i + 1],
							sba_problem->init_rot_+sba_problem->camera_index_[i]*4,
							sba_problem->init_rot_+nextcam_id*4,
							cameras+10,
							sba_problem->rolling_shutter_frame_time,
							sba_problem->rolling_shutter_capture_time,
							sba_problem->im_size[sba_problem->camera_index_[i]]));
					}else{
						// Bundler projection with horizontal rolling shutter
						cost_function =
							new AutoDiffCostFunction<ReprojectionErrorBundlerRadialRSVideoHorizontal, 2,9,9,3>(
							new ReprojectionErrorBundlerRadialRSVideoHorizontal(observations[2 * i + 0],
							observations[2 * i + 1],
							sba_problem->init_rot_+sba_problem->camera_index_[i]*4,
							sba_problem->init_rot_+nextcam_id*4,
							cameras+10,sba_problem->rolling_shutter_frame_time,
							sba_problem->rolling_shutter_capture_time,
							sba_problem->im_size[sba_problem->camera_index_[i]]));
					}
					break;
				case 51:
					num_camera_parameters = 9;
					if(sba_problem->rolling_shutter_orientation[sba_problem->camera_index_[i]]){
						// Bundler projection with vertical rolling shutter
						cost_function =
							new AutoDiffCostFunction<ReprojectionErrorBundlerRadialRSSingleWithCalib, 2,9,3,3,3>(
							new ReprojectionErrorBundlerRadialRSSingleWithCalib(observations+2 * i,
							sba_problem->init_rot_+sba_problem->camera_index_[i]*4,
							camera + 10,
							1
							));
					}else{
						num_camera_parameters = 9;
						// Bundler projection with horizontal rolling shutter
						cost_function =
							new AutoDiffCostFunction<ReprojectionErrorBundlerRadialRSSingleWithCalib, 2,9,3,3,3>(
							new ReprojectionErrorBundlerRadialRSSingleWithCalib(observations+2 * i,
							sba_problem->init_rot_ + sba_problem->camera_index_[i] * 4,
							camera + 10,
							0
							));
					}
					break;
				case 52:
					num_camera_parameters = 9;
					if (sba_problem->rolling_shutter_orientation[sba_problem->camera_index_[i]]){
						// Bundler projection with vertical rolling shutter
						cost_function =
							new AutoDiffCostFunction<ReprojectionErrorBundlerRadialRSSingleCayleyWithCalib, 2, 9, 3, 3, 3>(
							new ReprojectionErrorBundlerRadialRSSingleCayleyWithCalib(observations+2 * i ,
							sba_problem->init_rot_ + sba_problem->camera_index_[i] * 4,
							camera + 10,
							1
							));
					}
					else{
						num_camera_parameters = 9;
						// Bundler projection with horizontal rolling shutter
						cost_function =
							new AutoDiffCostFunction<ReprojectionErrorBundlerRadialRSSingleCayleyWithCalib, 2, 9, 3, 3, 3>(
							new ReprojectionErrorBundlerRadialRSSingleCayleyWithCalib(observations+2 * i,
							sba_problem->init_rot_ + sba_problem->camera_index_[i] * 4,
							camera + 10,
							0
							));
					}
					break;
				case 53:
					num_camera_parameters = 6;
					if (sba_problem->rolling_shutter_orientation[sba_problem->camera_index_[i]]){
						cost_function =
							new AutoDiffCostFunction<ReprojectionErrorBundlerRadialRSSingleShared, 2, 6, 3, 3, 3, 3>(
							new ReprojectionErrorBundlerRadialRSSingleShared(observations+2 * i,
							sba_problem->init_rot_ + sba_problem->camera_index_[i] * 4,
							camera + 10,
							1
							));
					}
					else{
						cost_function =
							new AutoDiffCostFunction<ReprojectionErrorBundlerRadialRSSingleShared, 2, 6, 3, 3, 3, 3>(
							new ReprojectionErrorBundlerRadialRSSingleShared(observations+2 * i,
							sba_problem->init_rot_ + sba_problem->camera_index_[i] * 4,
							camera + 10,
							0
							));
					}
					break;
				case 54:
					if (sba_problem->rolling_shutter_orientation[sba_problem->camera_index_[i]]){
						cost_function =
							new AutoDiffCostFunction<ReprojectionErrorBundlerRadialRSSingleCayleyShared, 2, 6, 3, 3, 3, 3>(
							new ReprojectionErrorBundlerRadialRSSingleCayleyShared(observations+2 * i,
							sba_problem->init_rot_ + sba_problem->camera_index_[i] * 4,
							camera + 10,
							1
							));
					}
					else{
						cost_function =
							new AutoDiffCostFunction<ReprojectionErrorBundlerRadialRSSingleCayleyShared, 2, 6, 3, 3, 3, 3>(
							new ReprojectionErrorBundlerRadialRSSingleCayleyShared(observations+2 * i,
							sba_problem->init_rot_ + sba_problem->camera_index_[i] * 4,
							camera + 10,
							0
							));
					}
					break;
				case 55:
					cost_function =
						new AutoDiffCostFunction<ReprojectionErrorBundlerRadialRSSingleUpVectorWithCalib, 2, 6, 3, 3, 1, 3>(
						new ReprojectionErrorBundlerRadialRSSingleUpVectorWithCalib(observations+2 * i,
						sba_problem->init_rot_ + sba_problem->camera_index_[i] * 4,
						camera + 10,
						1
						));
					break; 
				case 56:
					if (sba_problem->rolling_shutter_orientation[sba_problem->camera_index_[i]]){
						cost_function =
							new AutoDiffCostFunction<ReprojectionErrorBundlerRadialRSSingleUpVectorShared, 2, 3, 3, 3, 3, 1, 3>(
							new ReprojectionErrorBundlerRadialRSSingleUpVectorShared(observations+2 * i ,
							sba_problem->init_rot_ + sba_problem->camera_index_[i] * 4,
							camera + 10,
							1
							));
					}
					else{
						cost_function =
							new AutoDiffCostFunction<ReprojectionErrorBundlerRadialRSSingleUpVectorShared, 2, 3, 3, 3, 3, 1, 3>(
							new ReprojectionErrorBundlerRadialRSSingleUpVectorShared(observations+2 * i,
							sba_problem->init_rot_ + sba_problem->camera_index_[i] * 4,
							camera + 10,
							0
							));
					}
					break;
				case 57:
					if (sba_problem->rolling_shutter_orientation[sba_problem->camera_index_[i]]){
						cost_function =
							new AutoDiffCostFunction<ReprojectionErrorBundlerRadialRSSingleUpVectorVarShared, 2, 6, 3, 3, 3, 1, 3>(
							new ReprojectionErrorBundlerRadialRSSingleUpVectorVarShared(observations+2 * i,
							sba_problem->init_rot_ + sba_problem->camera_index_[i] * 4,
							camera + 10,
							1
							));
					}
					else{
						cost_function =
							new AutoDiffCostFunction<ReprojectionErrorBundlerRadialRSSingleUpVectorVarShared, 2, 6, 3, 3, 3, 1, 3>(
							new ReprojectionErrorBundlerRadialRSSingleUpVectorVarShared(observations+2 * i,
							sba_problem->init_rot_ + sba_problem->camera_index_[i] * 4,
							camera + 10,
							0
							));
					}
					break;
				case 58:
					if (sba_problem->rolling_shutter_orientation[sba_problem->camera_index_[i]]){
						cost_function =
							new AutoDiffCostFunction<ReprojectionErrorBundlerRadialRSSingleSharedWithoutWx, 2, 6, 3, 3, 2, 3>(
							new ReprojectionErrorBundlerRadialRSSingleSharedWithoutWx(observations + 2 * i,
							sba_problem->init_rot_ + sba_problem->camera_index_[i] * 4,
							camera + 10,
							*(sba_problem->rolling_shutter_w + sba_problem->camera_index_[i] * 3),
							1
							));
					}
					else{
						cost_function =
							new AutoDiffCostFunction<ReprojectionErrorBundlerRadialRSSingleSharedWithoutWx, 2, 6, 3, 3, 2, 3>(
							new ReprojectionErrorBundlerRadialRSSingleSharedWithoutWx(observations + 2 * i,
							sba_problem->init_rot_ + sba_problem->camera_index_[i] * 4,
							camera + 10,
							*(sba_problem->rolling_shutter_w + sba_problem->camera_index_[i] * 3),
							0
							));
					}
					break;
				case 59:
					if (sba_problem->rolling_shutter_orientation[sba_problem->camera_index_[i]]){
						// Bundler projection with vertical rolling shutter
						cost_function =
							new AutoDiffCostFunction<ReprojectionErrorRSForssen, 2, 9, 3, 3, 3>(
							new ReprojectionErrorRSForssen(observations + 2 * i,
							sba_problem->init_rot_ + sba_problem->camera_index_[i] * 4,
							camera + 10,
							1
							));
					}
					else{
						// Bundler projection with horizontal rolling shutter
						cost_function =
							new AutoDiffCostFunction<ReprojectionErrorRSForssen, 2, 9, 3, 3, 3>(
							new ReprojectionErrorRSForssen(observations + 2 * i,
							sba_problem->init_rot_ + sba_problem->camera_index_[i] * 4,
							camera + 10, 
							0
							));
					}
					break;

				case 60:
					num_camera_parameters = 9;
					if (sba_problem->rolling_shutter_orientation[sba_problem->camera_index_[i]]){
						// Bundler projection with vertical rolling shutter
						cost_function =
							new AutoDiffCostFunction<ReprojectionErrorBundlerRadialRSSingleEAXWithCalib, 2, 9, 3, 3, 3>(
							new ReprojectionErrorBundlerRadialRSSingleEAXWithCalib(observations + 2 * i,
							sba_problem->init_rot_ + sba_problem->camera_index_[i] * 4,
							camera + 10,
							1
							));
					}
					else{
						num_camera_parameters = 9;
						// Bundler projection with horizontal rolling shutter
						cost_function =
							new AutoDiffCostFunction<ReprojectionErrorBundlerRadialRSSingleEAXWithCalib, 2, 9, 3, 3, 3>(
							new ReprojectionErrorBundlerRadialRSSingleEAXWithCalib(observations + 2 * i,
							sba_problem->init_rot_ + sba_problem->camera_index_[i] * 4,
							camera + 10,
							0
							));
					}
					break;
					case 100:
						// Bundler projection with vertical rolling shutter
						cost_function =
							new AutoDiffCostFunction<KalibrRadtanReprojectionError, 2, 21, 3>(
							new KalibrRadtanReprojectionError(observations[2 * i], 
							observations[2 * i +1]
							));
					break;
				}





				// If enabled use Huber's loss function.
				LossFunction* loss_function = sba_problem->robustify ? new HuberLoss(1.0) : NULL;

				// Each observation correponds to a pair of a camera and a point
				// which are identified by camera_index_[i] and point_index_[i]
				// respectively.

				double* point = points + point_block_size * sba_problem->point_index_[i];
				if (sba_problem->covariance_mask != NULL){
					if (sba_problem->covariance_mask[sba_problem->num_cameras_ +sba_problem->point_index_[i]] > 0){
						sba_problem->covariance_blocks.push_back(make_pair(point, point));
						sba_problem->covariance_mask[sba_problem->num_cameras_ + sba_problem->point_index_[i]]--;
					}
				}
				double *q, *t, *fx, *ar, *x0y0, *s, *r;
				switch(sba_problem->proj_func_){
				case 0:
					
					q = camera + 1;
					t = camera + 4;
					fx = camera + 9;
					ar = camera + 12;
					x0y0 = camera + 10;
					s = camera + 13;
					r = camera + 7;
					if (sba_problem->shared_parmask != NULL){
						if (sba_problem->shared_parmask[sba_problem->camera_index_[i] * 7]>=0){
							//shared q
							q = cameras + sba_problem->shared_parmask[sba_problem->camera_index_[i] * 7] * sba_problem->num_camera_parameters_ + 1;
						}
						
						if (sba_problem->shared_parmask[sba_problem->camera_index_[i] * 7 + 1] >= 0){
							//shared t
							t = cameras + sba_problem->shared_parmask[sba_problem->camera_index_[i] * 7+1] * sba_problem->num_camera_parameters_ + 4;
						}
						
						if (sba_problem->shared_parmask[sba_problem->camera_index_[i] * 7 + 2] >= 0){
							//shared fx
							fx = cameras + sba_problem->shared_parmask[sba_problem->camera_index_[i] * 7+2] * sba_problem->num_camera_parameters_ + 9;
						}
						
						if (sba_problem->shared_parmask[sba_problem->camera_index_[i] * 7 + 3] >= 0){
							//shared ar
							ar = cameras + sba_problem->shared_parmask[sba_problem->camera_index_[i] * 7+3] * sba_problem->num_camera_parameters_ + 12;
						}
						
						if (sba_problem->shared_parmask[sba_problem->camera_index_[i] * 7 + 4] >= 0){
							//shared x0y0
							x0y0 = cameras + sba_problem->shared_parmask[sba_problem->camera_index_[i] * 7+4] * sba_problem->num_camera_parameters_ + 10;
						}
						
						if (sba_problem->shared_parmask[sba_problem->camera_index_[i] * 7 + 5] >= 0){
							//shared s
							s = cameras + sba_problem->shared_parmask[sba_problem->camera_index_[i] * 7+5] * sba_problem->num_camera_parameters_ + 13;
						}

						if (sba_problem->shared_parmask[sba_problem->camera_index_[i] * 7 + 6] >= 0){
							//shared r
							r = cameras + sba_problem->shared_parmask[sba_problem->camera_index_[i] * 7+6] * sba_problem->num_camera_parameters_ + 7;
						}
						
					}
					problem->AddResidualBlock(cost_function,
						loss_function,
						q,
						t,
						fx,
						x0y0,
						ar,
						s,
						r,
						point);
					if (sba_problem->covariance_mask != NULL){
						if (sba_problem->covariance_mask[sba_problem->camera_index_[i]] > 0){
							sba_problem->covariance_blocks.push_back(make_pair(q, q));
							sba_problem->covariance_blocks.push_back(make_pair(t, t));
							sba_problem->covariance_blocks.push_back(make_pair(fx, fx));
							sba_problem->covariance_blocks.push_back(make_pair(x0y0, x0y0));
							sba_problem->covariance_blocks.push_back(make_pair(ar, ar));
							sba_problem->covariance_blocks.push_back(make_pair(s, s));
							sba_problem->covariance_blocks.push_back(make_pair(r, r));
							sba_problem->covariance_mask[sba_problem->camera_index_[i]]--;
						}
					}
					break;
				case 1:
				case 2:
				case 3:
				case 4:
				case 5:
				case 6:
				case 7:
				case 8:
				case 9:
				case 10:
				case 11:
				case 20:
				case 21:
					problem->AddResidualBlock(cost_function,
						loss_function,
						camera+1,
						point);
					if (sba_problem->covariance_mask != NULL){
						if (sba_problem->covariance_mask[sba_problem->camera_index_[i]] > 0){
							sba_problem->covariance_blocks.push_back(make_pair(camera + 1, camera + 1));
							sba_problem->covariance_mask[sba_problem->camera_index_[i]]--;
						}
					}
					break;
				case 12:
					problem->AddResidualBlock(cost_function,
						loss_function,
						camera,
						point);
					break;
				case 22:
				case 23:
					problem->AddResidualBlock(cost_function,
						loss_function,
						camera,
						camera+4,
						point);
					if (sba_problem->covariance_mask != NULL){
						if (sba_problem->covariance_mask[sba_problem->camera_index_[i]] > 0){
							sba_problem->covariance_blocks.push_back(make_pair(camera + 1, camera + 1));
							sba_problem->covariance_mask[sba_problem->camera_index_[i]]--;
						}
					}
					break;
				case 24:
				case 25:
					problem->AddResidualBlock(cost_function,
						loss_function,
						camera,
						camera+4,
						point);
					problem->SetParameterization(camera,quaternion_parameterization);
					if (sba_problem->covariance_mask != NULL){
						if (sba_problem->covariance_mask[sba_problem->camera_index_[i]] > 0){
							sba_problem->covariance_blocks.push_back(make_pair(camera, camera));
							sba_problem->covariance_blocks.push_back(make_pair(camera + 4, camera + 4));
							sba_problem->covariance_mask[sba_problem->camera_index_[i]]--;
						}
					}
					break;
				case 30:
				case 31:
					problem->AddResidualBlock(cost_function,
						loss_function,
						camera,
						point);
					if (sba_problem->covariance_mask != NULL){
						if (sba_problem->covariance_mask[sba_problem->camera_index_[i]] > 0){
							sba_problem->covariance_blocks.push_back(make_pair(camera, camera));
							sba_problem->covariance_mask[sba_problem->camera_index_[i]]--;
						}
					}
				case 40:
					problem->AddResidualBlock(cost_function,
						loss_function,
						camera+1,
						point);
					if (sba_problem->covariance_mask != NULL){
						if (sba_problem->covariance_mask[sba_problem->camera_index_[i]] > 0){
							sba_problem->covariance_blocks.push_back(make_pair(camera + 1, camera + 1));
							sba_problem->covariance_mask[sba_problem->camera_index_[i]]--;
						}
					}
					break;
				case 41:
					problem->AddResidualBlock(cost_function,
						loss_function,
						camera+1, 
						cameras+7,
						point);
					if (sba_problem->covariance_mask != NULL){
						if (sba_problem->covariance_mask[sba_problem->camera_index_[i]] > 0){
							sba_problem->covariance_blocks.push_back(make_pair(camera + 1, camera + 1));
							sba_problem->covariance_blocks.push_back(make_pair(camera + 7, camera + 7));
							sba_problem->covariance_mask[sba_problem->camera_index_[i]]--;
						}
					}
					break;
				case 42:
					problem->AddResidualBlock(cost_function,
						loss_function,
						camera + 1,
						camera + 7,
						point);
						problem->SetParameterBlockConstant(camera + 7);
						break;
				case 50:

					problem->AddResidualBlock(cost_function,
						loss_function,
						camera+1, 
						next_camera+1, //next camera
						point);
					if (sba_problem->covariance_mask != NULL){
						if (sba_problem->covariance_mask[sba_problem->camera_index_[i]] > 0){
							sba_problem->covariance_blocks.push_back(make_pair(camera + 1, camera + 1));
							sba_problem->covariance_mask[sba_problem->camera_index_[i]]--;
						}
					}
					break;
				case 51:
				case 52:
				case 59:
				case 60:
					problem->AddResidualBlock(cost_function,
						loss_function,
						camera + 1,
						sba_problem->rolling_shutter_v +sba_problem->camera_index_[i] * 3,
						sba_problem->rolling_shutter_w + sba_problem->camera_index_[i] * 3,
						point);
					if (sba_problem->covariance_mask != NULL){
						if (sba_problem->covariance_mask[sba_problem->camera_index_[i]] > 0){
							sba_problem->covariance_blocks.push_back(make_pair(camera + 1, camera + 1));
							sba_problem->covariance_mask[sba_problem->camera_index_[i]]--;
						}
					}
					
					break;
				case 53:
					problem->AddResidualBlock(cost_function,
						loss_function,
						camera + 1,
						camera + 7,
						sba_problem->rolling_shutter_v + sba_problem->camera_index_[i] * 3,
						sba_problem->rolling_shutter_w + sba_problem->camera_index_[i] * 3,
						point);
					problem->SetParameterBlockConstant(camera+7);
					break;
				case 54:
					problem->AddResidualBlock(cost_function,
						loss_function,
						camera + 1,
						camera + 7,
						sba_problem->rolling_shutter_v + sba_problem->camera_index_[i] * 3,
						sba_problem->rolling_shutter_w + sba_problem->camera_index_[i] * 3,
						point);
					problem->SetParameterBlockConstant(camera + 7);
					break;
				case 55:
					problem->AddResidualBlock(cost_function,
						loss_function,
						camera + 4,
						sba_problem->rolling_shutter_v + sba_problem->camera_index_[i] * 3,
						sba_problem->rolling_shutter_w + sba_problem->camera_index_[i] * 3,
						sba_problem->q + +sba_problem->camera_index_[i],
						point);
					break;
				case 56:
					problem->AddResidualBlock(cost_function,
						loss_function,
						camera + 4,
						camera + 7,
						sba_problem->rolling_shutter_v + sba_problem->camera_index_[i] * 3,
						sba_problem->rolling_shutter_w + sba_problem->camera_index_[i] * 3,
						sba_problem->q + +sba_problem->camera_index_[i],
						point);
					problem->SetParameterBlockConstant(camera + 7);
					break;
				case 57:
					problem->AddResidualBlock(cost_function,
						loss_function,
						camera + 1,
						camera + 7,
						sba_problem->rolling_shutter_v + sba_problem->camera_index_[i] * 3,
						sba_problem->rolling_shutter_w + sba_problem->camera_index_[i] * 3,
						sba_problem->q + +sba_problem->camera_index_[i],
						point);
					problem->SetParameterBlockConstant(camera + 7);
					break;
				case 58:
					problem->AddResidualBlock(cost_function,
						loss_function,
						camera + 1,
						camera + 7,
						sba_problem->rolling_shutter_v + sba_problem->camera_index_[i] * 3,
						sba_problem->rolling_shutter_w + sba_problem->camera_index_[i] * 3+1,
						sba_problem->q + +sba_problem->camera_index_[i],
						point);
					problem->SetParameterBlockConstant(camera + 7);
					break;
				case 100:
					problem->AddResidualBlock(cost_function,
						loss_function,
						camera,
						point);
					break;
				}

				if(sba_problem->proj_func_==0){
					if (sba_problem->fixed_parmask != NULL){
						if (sba_problem->fixed_parmask[sba_problem->camera_index_[i] * 7])problem->SetParameterBlockConstant(q); //constant q
						if (sba_problem->fixed_parmask[sba_problem->camera_index_[i] * 7 + 1])problem->SetParameterBlockConstant(t); //constant t
						if (sba_problem->fixed_parmask[sba_problem->camera_index_[i] * 7 + 2])problem->SetParameterBlockConstant(fx); //constant fx
						if (sba_problem->fixed_parmask[sba_problem->camera_index_[i] * 7 + 3])problem->SetParameterBlockConstant(ar); //constant ar
						if (sba_problem->fixed_parmask[sba_problem->camera_index_[i] * 7 + 4])problem->SetParameterBlockConstant(x0y0); //constant x0y0
						if (sba_problem->fixed_parmask[sba_problem->camera_index_[i] * 7 + 5])problem->SetParameterBlockConstant(s); //constant s
						if (sba_problem->fixed_parmask[sba_problem->camera_index_[i] * 7 + 6])problem->SetParameterBlockConstant(r); //constant r
					}
				}

				if (sba_problem->constant_cameras[sba_problem->camera_index_[i]])problem->SetParameterBlockConstant(camera + 1);
				if (sba_problem->constant_points[sba_problem->point_index_[i]])problem->SetParameterBlockConstant(point);





			}

			if (sba_problem->fixed_parmask != NULL && sba_problem->proj_func_!=0){
				for (int j = 0; j < sba_problem->num_cameras_; j++){
					ceres::SubsetParameterization *constant_parameterization = NULL;
					double * camera = sba_problem->parameters_ + camera_block_size*j+1;
					std::vector<int> constant_pars;
					for (int i = 0; i < 9; i++)
					{
						if (sba_problem->fixed_parmask[j * 17 + i]){
							constant_pars.push_back(i);
						}
					}
					constant_parameterization = new ceres::SubsetParameterization(9, constant_pars);
					problem->SetParameterization(camera, constant_parameterization);
				}
			}

			if(sba_problem->constrained_cameras){
				CostFunction* cost_function;
				for(int i=0;i<sba_problem->num_cameras_;i++){
					double * camera = sba_problem->parameters_+camera_block_size*i;
					double * weights = sba_problem->camera_constraint_weights+camera_block_size * i;
					double * constraints = sba_problem->camera_constraints+camera_block_size * i;

					switch(sba_problem->proj_func_){
					case 0:
						//q
						cost_function =
							new AutoDiffCostFunction<CameraConstraintsError, 3, 3>(
							new CameraConstraintsError(constraints + 1, weights+1, 3));
						problem->AddResidualBlock(cost_function,
							NULL,
							camera + 1);
						//t
						cost_function =
							new AutoDiffCostFunction<CameraConstraintsError, 3, 3>(
							new CameraConstraintsError(constraints + 4, weights + 4, 3));
						problem->AddResidualBlock(cost_function,
							NULL,
							camera + 4);
						//r
						cost_function =
							new AutoDiffCostFunction<CameraConstraintsError, 2, 2>(
							new CameraConstraintsError(constraints + 7, weights + 7, 2));
						problem->AddResidualBlock(cost_function,
							NULL,
							camera + 7);
						//fx
						cost_function =
							new AutoDiffCostFunction<CameraConstraintsError, 1, 1>(
							new CameraConstraintsError(constraints + 9, weights + 9, 1));
						problem->AddResidualBlock(cost_function,
							NULL,
							camera + 9);
						//x0y0
						cost_function =
							new AutoDiffCostFunction<CameraConstraintsError, 2, 2>(
							new CameraConstraintsError(constraints + 10, weights + 10, 2));
						problem->AddResidualBlock(cost_function,
							NULL,
							camera + 10);
						//ar
						cost_function =
							new AutoDiffCostFunction<CameraConstraintsError, 1, 1>(
							new CameraConstraintsError(constraints + 12, weights + 12, 1));
						problem->AddResidualBlock(cost_function,
							NULL,
							camera + 12);
						//s
						cost_function =
							new AutoDiffCostFunction<CameraConstraintsError, 1, 1>(
							new CameraConstraintsError(constraints + 13, weights + 13, 1));
						problem->AddResidualBlock(cost_function,
							NULL,
							camera + 13);
						break;
					case 1:
						cost_function =
							new AutoDiffCostFunction<CameraConstraintsError,11,11>(
							new CameraConstraintsError(constraints+1,weights,11));
						break;
					case 2:
						cost_function =
							new AutoDiffCostFunction<CameraConstraintsError,6,6>(
							new CameraConstraintsError(constraints + 1, weights+1, 6));
						break;
					case 3:
						cost_function =
							new AutoDiffCostFunction<CameraConstraintsError,9,9>(
							new CameraConstraintsError(constraints + 1, weights + 1, 9));
						break;
					case 4:
						cost_function =
							new AutoDiffCostFunction<CameraConstraintsError,10,10>(
							new CameraConstraintsError(constraints + 1, weights + 1, 10));
						break;
					case 5:
						cost_function =
							new AutoDiffCostFunction<CameraConstraintsError,7,7>(
							new CameraConstraintsError(constraints + 1, weights + 1, 7));
						break;
					case 6:
						cost_function =
							new AutoDiffCostFunction<CameraConstraintsError,13,13>(
							new CameraConstraintsError(constraints + 1, weights + 1, 13));
						break;
					case 7:
						cost_function =
							new AutoDiffCostFunction<CameraConstraintsError,6,6>(
							new CameraConstraintsError(constraints + 1, weights + 1, 6));
						break;
					case 8:
						cost_function =
							new AutoDiffCostFunction<CameraConstraintsError,11,11>(
							new CameraConstraintsError(constraints + 1, weights + 1, 11));
						break;
					case 9:
						cost_function =
							new AutoDiffCostFunction<CameraConstraintsError,12,12>(
							new CameraConstraintsError(constraints + 1, weights + 1, 12));
						break;
					case 10:
						cost_function =
							new AutoDiffCostFunction<CameraConstraintsError,9,9>(
							new CameraConstraintsError(constraints + 1, weights + 1, 9));
						break;
					case 11:
						cost_function =
							new AutoDiffCostFunction<CameraConstraintsError, 9, 9>(
							new CameraConstraintsError(constraints + 1, weights + 1, 9));
						break;
					case 12:
						cost_function =
							new AutoDiffCostFunction<CameraConstraintsError, 12, 12>(
							new CameraConstraintsError(constraints, weights, 12));
						break;
					case 20:
						cost_function =
							new AutoDiffCostFunction<CameraConstraintsError,7,7>(
							new CameraConstraintsError(constraints + 1, weights + 1, 7));
						break;
					case 21:
						cost_function =
							new AutoDiffCostFunction<CameraConstraintsError,9,9>(
							new CameraConstraintsError(constraints + 1, weights + 1, 9));
						break;
					case 22:
					case 24:
						cost_function =
							new AutoDiffCostFunction<CameraConstraintsError,4,4>(
							new CameraConstraintsError(constraints + 4, weights + 4, 4));
						break;
					case 23:
					case 25:
						cost_function =
							new AutoDiffCostFunction<CameraConstraintsError,6,6>(
							new CameraConstraintsError(constraints + 4, weights + 4, 6));
						break;
					case 30: 
						cost_function =
							new AutoDiffCostFunction<CameraConstraintsError,7,7>(
							new CameraConstraintsError(constraints, weights, 7));
						break;
					case 31: 
						cost_function =
							new AutoDiffCostFunction<CameraConstraintsError,9,9>(
							new CameraConstraintsError(constraints, weights, 9));
						break;
					case 40:
						cost_function =
							new AutoDiffCostFunction<CameraConstraintsError,9,9>(
							new CameraConstraintsError(constraints+1, weights+1 , 9));
						break;
					case 41:
					case 42:
						cost_function =
							new AutoDiffCostFunction<CameraConstraintsError,3,3>(
							new CameraConstraintsError(constraints+7,weights+7,3));
						break;
					case 50:
						cost_function =
							new AutoDiffCostFunction<CameraConstraintsError,9,9>(
							new CameraConstraintsError(constraints + 1, weights + 1, 9));
						break;
					case 51:
					case 52:					
					case 59:
					case 60:
						cost_function =
							new AutoDiffCostFunction<CameraConstraintsError, 9, 9>(
							new CameraConstraintsError(constraints + 1, weights + 1, 9));
						problem->AddResidualBlock(cost_function,
							NULL,
							camera + 1);
						if (sba_problem->rolling_shutter_v_constraints != NULL){
							cost_function = new AutoDiffCostFunction<CameraConstraintsError, 3, 3>(
								new CameraConstraintsError(sba_problem->rolling_shutter_v_constraints + i * 3, sba_problem->rolling_shutter_v_weights + i * 3, 3));
							problem->AddResidualBlock(cost_function,
								NULL,
								sba_problem->rolling_shutter_v + i * 3);
							cost_function = new AutoDiffCostFunction<CameraConstraintsError, 3, 3>(
								new CameraConstraintsError(sba_problem->rolling_shutter_w_constraints + i * 3, sba_problem->rolling_shutter_w_weights + i * 3, 3));
							problem->AddResidualBlock(cost_function,
								NULL,
								sba_problem->rolling_shutter_w + i * 3);

						}
						break;
					case 53:
					case 54:
						if (sba_problem->rolling_shutter_v_constraints != NULL){
							cost_function = new AutoDiffCostFunction<CameraConstraintsError, 3, 3>(
								new CameraConstraintsError(sba_problem->rolling_shutter_v_constraints + i * 3, sba_problem->rolling_shutter_v_weights + i * 3, 3));
							problem->AddResidualBlock(cost_function,
								NULL,
								sba_problem->rolling_shutter_v + i * 3);
							cost_function = new AutoDiffCostFunction<CameraConstraintsError, 3, 3>(
								new CameraConstraintsError(sba_problem->rolling_shutter_w_constraints + i * 3, sba_problem->rolling_shutter_w_weights + i * 3, 3));
							problem->AddResidualBlock(cost_function,
								NULL,
								sba_problem->rolling_shutter_w + i * 3);

						}
						break;
					case 55:
						cost_function =
							new AutoDiffCostFunction<CameraConstraintsError, 6, 6>(
							new CameraConstraintsError(constraints + 4, weights + 4, 6));
						problem->AddResidualBlock(cost_function,
							NULL,
							camera + 4);
						if (sba_problem->rolling_shutter_v_constraints != NULL){
							cost_function = new AutoDiffCostFunction<CameraConstraintsError, 3, 3>(
								new CameraConstraintsError(sba_problem->rolling_shutter_v_constraints + i * 3, sba_problem->rolling_shutter_v_weights + i * 3, 3));
							problem->AddResidualBlock(cost_function,
								NULL,
								sba_problem->rolling_shutter_v + i * 3);
							cost_function = new AutoDiffCostFunction<CameraConstraintsError, 3, 3>(
								new CameraConstraintsError(sba_problem->rolling_shutter_w_constraints + i * 3, sba_problem->rolling_shutter_w_weights + i * 3, 3));
							problem->AddResidualBlock(cost_function,
								NULL,
								sba_problem->rolling_shutter_w + i * 3);

						}
						break;
					case 56:
						if (sba_problem->rolling_shutter_v_constraints != NULL){
							cost_function = new AutoDiffCostFunction<CameraConstraintsError, 3, 3>(
								new CameraConstraintsError(sba_problem->rolling_shutter_v_constraints + i * 3, sba_problem->rolling_shutter_v_weights + i * 3, 3));
							problem->AddResidualBlock(cost_function,
								NULL,
								sba_problem->rolling_shutter_v + i * 3);
							cost_function = new AutoDiffCostFunction<CameraConstraintsError, 3, 3>(
								new CameraConstraintsError(sba_problem->rolling_shutter_w_constraints + i * 3, sba_problem->rolling_shutter_w_weights + i * 3, 3));
							problem->AddResidualBlock(cost_function,
								NULL,
								sba_problem->rolling_shutter_w + i * 3);

						}
						break;
					case 57:
						if (sba_problem->rolling_shutter_v_constraints != NULL){
							cost_function = new AutoDiffCostFunction<CameraConstraintsError, 3, 3>(
								new CameraConstraintsError(sba_problem->rolling_shutter_v_constraints + i * 3, sba_problem->rolling_shutter_v_weights + i * 3, 3));
							problem->AddResidualBlock(cost_function,
								NULL,
								sba_problem->rolling_shutter_v + i * 3);
							cost_function = new AutoDiffCostFunction<CameraConstraintsError, 3, 3>(
								new CameraConstraintsError(sba_problem->rolling_shutter_w_constraints + i * 3, sba_problem->rolling_shutter_w_weights + i * 3, 3));
							problem->AddResidualBlock(cost_function,
								NULL,
								sba_problem->rolling_shutter_w + i * 3);

						}
						break;
					case 58:
						if (sba_problem->rolling_shutter_v_constraints != NULL){
							cost_function = new AutoDiffCostFunction<CameraConstraintsError, 3, 3>(
								new CameraConstraintsError(sba_problem->rolling_shutter_v_constraints + i * 3, sba_problem->rolling_shutter_v_weights + i * 3, 3));
							problem->AddResidualBlock(cost_function,
								NULL,
								sba_problem->rolling_shutter_v + i * 3);
							cost_function = new AutoDiffCostFunction<CameraConstraintsError, 2, 2>(
								new CameraConstraintsError(sba_problem->rolling_shutter_w_constraints + i * 3+1, sba_problem->rolling_shutter_w_weights + i * 3+1, 2));
							problem->AddResidualBlock(cost_function,
								NULL,
								sba_problem->rolling_shutter_w + i * 3+1);

						}
						break;

					}

					switch(sba_problem->proj_func_){
					case 0:
						break;
					case 1:
					case 2:
					case 3:
					case 4:
					case 5:
					case 6:
					case 7:
					case 8:
					case 9:
					case 10:
					case 11:
					case 20:
					case 21:
						problem->AddResidualBlock(cost_function,
							NULL,
							camera+1);
						break;
					case 12:
						problem->AddResidualBlock(cost_function,
							NULL,
							camera);
						break;
					case 22:
					case 23:
						problem->AddResidualBlock(cost_function,
							NULL,
							camera+4);
						break;
					case 24:
					case 25:
						problem->AddResidualBlock(cost_function,
							NULL,
							camera+4);
						break;
					case 30:
					case 31:
						problem->AddResidualBlock(cost_function,
							NULL,
							camera);
					case 40:
						problem->AddResidualBlock(cost_function,
							NULL,
							camera+1);
						break;
					case 41:
					case 42:
						problem->AddResidualBlock(cost_function,
							NULL,
							cameras+7);
						break;
					case 50:
						problem->AddResidualBlock(cost_function,
							NULL,
							camera+1);
						break;
					case 51:
					case 52:
					case 60:
												break;
					}

				}
			}

			if (sba_problem->n_ref_pts_ > 0){
				addReferencePoints(sba_problem, problem);
			}

			if (sba_problem->n_ref_cams_ > 0){
				addReferenceCameras(sba_problem, problem);
			}

		}



		//void BuildPancamProblem(BAProblem* sba_problem, Problem* problem) {
		//  const int point_block_size = sba_problem->num_point_parameters_;
		//  const int camera_block_size = sba_problem->num_camera_parameters_;
		//  double* points = sba_problem->mutable_points();
		//  double* cameras = sba_problem->mutable_cameras();
		//  int num_calibration_parameters;
		//
		//    if(sba_problem->estimate_radial_){
		//  	  num_calibration_parameters = 7;
		//    }else{
		//  	  num_calibration_parameters = 5;
		//  }
		//
		//  // Observations is 2*num_observations long array observations =
		//  // [u_1, u_2, ... , u_n], where each u_i is two dimensional, the x
		//  // and y positions of the observation.
		//  const double* observations = sba_problem->observations_;
		//  //save the initial rotations
		//  double * initrots = new double[sba_problem->num_cameras_*8];
		//  for(int i=0;i<sba_problem->num_cameras_;i++){
		//	  double * q = cameras+i*camera_block_size + num_calibration_parameters;
		//	  double * qc = cameras+i*camera_block_size + num_calibration_parameters+8;
		//	  initrots[i*8]=q[0];
		//	  initrots[i*8+1]=q[1];
		//	  initrots[i*8+2]=q[2];
		//	  initrots[i*8+3]=q[3];
		//	  //initial rotation for camera circle
		//	  initrots[i*8+4]=qc[0];
		//	  initrots[i*8+5]=qc[1];
		//	  initrots[i*8+6]=qc[2];
		//	  initrots[i*8+7]=qc[3];
		//	  //set local rotations to zero
		//	  q[0]=1;
		//	  q[1]=0;
		//	  q[2]=0;
		//	  q[3]=0;
		//	  qc[0]=1;
		//	  qc[1]=0;
		//	  qc[2]=0;
		//	  qc[3]=0;
		//
		//
		//  }
		//
		//  for (int i = 0; i < sba_problem->num_observations_; ++i) {
		//    CostFunction* cost_function;
		//    // Each Residual block takes a point and a camera as input and
		//    // outputs a 2 dimensional residual.
		//    double* camera = cameras + camera_block_size * sba_problem->camera_index_[i];
		//    if(sba_problem->estimate_radial_){
		//    	cost_function =
		//			  new AutoDiffCostFunction<PancamReprojectionErrorKQTRadial, 2, 7,4,3,3,1,3>(
		//				  new PancamReprojectionErrorKQTRadial(observations[2 * i + 0],
		//											   observations[2 * i + 1],
		//											   initrots+sba_problem->camera_index_[i]*8,
		//											   initrots+sba_problem->camera_index_[i]*8+4,
		//											   camera[3],camera[4]));
		//    }else{
		//    	cost_function =
		//			  new AutoDiffCostFunction<PancamReprojectionErrorKQT, 2, 5,4,3,3,1,3>(
		//				  new PancamReprojectionErrorKQT(observations[2 * i + 0],
		//											   observations[2 * i + 1],
		//											   initrots+sba_problem->camera_index_[i]*8,
		//											   initrots+sba_problem->camera_index_[i]*8+4,
		//											   camera[3],camera[4]));
		//
		//    }
		//
		//    // If enabled use Huber's loss function.
		//    LossFunction* loss_function = FLAGS_robustify ? new HuberLoss(1.0) : NULL;
		//
		//    // Each observation correponds to a pair of a camera and a point
		//    // which are identified by camera_index_[i] and point_index_[i]
		//    // respectively.
		//
		//    double* point = points + point_block_size * sba_problem->point_index_[i];
		//
		//    double * Cc,*Qc,*rho;
		//    //shared parameters of the circle center
		//    int block_cluster = sba_problem->parmask_.block_cam_mask[0][sba_problem->camera_index_[i]];
		//    Cc = sba_problem->parmask_.block_pointers[0][block_cluster-1][0]+num_calibration_parameters+5;
		//    //shared parameters of the circle orientation
		//    block_cluster = sba_problem->parmask_.block_cam_mask[1][sba_problem->camera_index_[i]];
		//    Qc = sba_problem->parmask_.block_pointers[1][block_cluster-1][0]+num_calibration_parameters+9;
		//    //shared parameters of the circle radius
		//    block_cluster = sba_problem->parmask_.block_cam_mask[2][sba_problem->camera_index_[i]];
		//    rho = sba_problem->parmask_.block_pointers[2][block_cluster-1][0]+num_calibration_parameters+12;
		//
		//    problem->AddResidualBlock(cost_function,
		//    		 	 	 	 	  loss_function,
		//    		 	 	 	 	  camera,
		//    		 	 	 	 	  camera+num_calibration_parameters+1,
		//    		 	 	 	 	  Cc,
		//    		 	 	 	 	  Qc,
		//    		 	 	 	 	  rho,
		//    		 	 	 	 	  point);
		//
		//  }
		//
		//  delete [] initrots;
		//
		//}

		//void BuildMulticamProblem(BAProblem* sba_problem, Problem* problem) {
		//  const int point_block_size = sba_problem->num_point_parameters_;
		//  const int camera_block_size = sba_problem->num_camera_parameters_;
		//  double* points = sba_problem->mutable_points();
		//  double* cameras = sba_problem->mutable_cameras();
		//
		//  // Observations is 2*num_observations long array observations =
		//  // [u_1, u_2, ... , u_n], where each u_i is two dimensional, the x
		//  // and y positions of the observation.
		//  const double* observations = sba_problem->observations_;
		//  //save the initial rotations
		//  double * initrots = new double[sba_problem->num_cameras_*4];
		//
		//  int camid,clusterid,cam_clusterid;
		//  for (int i = 0; i < sba_problem->num_observations_; ++i) {
		//    CostFunction* cost_function;
		//
		//	camid = sba_problem->camera_index_[i];
		//	clusterid = sba_problem->cam_clusters_[camid];
		//	cam_clusterid = sba_problem->camids_[camid];
		//    double* camera = cameras + camera_block_size * camid;
		//    if(sba_problem->estimate_radial_){
		//		switch(sba_problem->proj_func_){
		//			case 200:
		//				cost_function =
		//				new AutoDiffCostFunction<ReprojectionErrorMulticamKQT, 2, 5,3,3,3,3 ,3>(
		//				new ReprojectionErrorMulticamKQT(observations[2 * i + 0],observations[2 * i + 1],sba_problem->init_rot_+cam_clusterid*4,sba_problem->init_rot_cluster_+clusterid*4));
		//				break;
		//		}
		//    }else{
		//		switch(sba_problem->proj_func_){
		//		case 200:
		//    		cost_function =
		//    		new AutoDiffCostFunction<ReprojectionErrorMulticamKQT, 2, 5,3,3,3,3 ,3>(
		//			new ReprojectionErrorMulticamKQT(observations[2 * i + 0], observations[2 * i + 1],sba_problem->init_rot_+cam_clusterid*4,sba_problem->init_rot_cluster_+clusterid*4));
		//			break;
		//		}
		//    }
		//
		//
		//    // If enabled use Huber's loss function.
		//    LossFunction* loss_function = FLAGS_robustify ? new HuberLoss(1.0) : NULL;
		//
		//    // Each observation correponds to a pair of a camera and a point
		//    // which are identified by camera_index_[i] and point_index_[i]
		//    // respectively.
		//
		//    double* point = points + point_block_size * sba_problem->point_index_[i];
		//
		//	
		//	problem->AddResidualBlock(cost_function,
		//    		 	 	 	 		loss_function,
		//    		 	 	 	 		camera+7,
		//								sba_problem->Rc_+cam_clusterid*3,
		//								sba_problem->tc_+cam_clusterid*3,
		//								sba_problem->Rs_+clusterid*3,
		//								sba_problem->ts_+clusterid*3,
		//    		 	 	 	 		point);
		//		
		//	
		//	//if(FLAGS_constant_K)problem->SetParameterBlockConstant(camera);
		//	/*if(FLAGS_constant_R)problem->SetParameterBlockConstant(q);*/
		//	if(FLAGS_constant_cameras)problem->SetParameterBlockConstant(camera+1);
		//	if(FLAGS_constant_points)problem->SetParameterBlockConstant(point);	
		//
		//
		//  }
		//  delete [] initrots;
		//
		//}


		//void BuildCalibProblem(BAProblem* sba_problem, Problem* problem) {
		//  const int point_block_size = sba_problem->num_point_parameters_;
		//  const int camera_block_size = sba_problem->num_camera_parameters_;
		//  double* points = sba_problem->mutable_points();
		//  double* cameras = sba_problem->mutable_cameras();
		//
		//  // Observations is 2*num_observations long array observations =
		//  // [u_1, u_2, ... , u_n], where each u_i is two dimensional, the x
		//  // and y positions of the observation.
		//  const double* observations = sba_problem->observations_;
		//  //save the initial rotations
		//  double * initrots = new double[sba_problem->num_cameras()*4];
		//  for(int i=0;i<sba_problem->num_cameras();i++){
		//	  double * q = cameras+i*camera_block_size;
		//	  initrots[i*4]=q[0];
		//	  initrots[i*4+1]=q[1];
		//	  initrots[i*4+2]=q[2];
		//	  initrots[i*4+3]=q[3];
		//	  sba_problem->init_rot_[i*4]=q[0];
		//	  sba_problem->init_rot_[i*4+1]=q[1];
		//	  sba_problem->init_rot_[i*4+2]=q[2];
		//	  sba_problem->init_rot_[i*4+3]=q[3];
		//	  //set local rotations to zero
		//	  q[0]=1;
		//	  q[1]=0;
		//	  q[2]=0;
		//	  q[3]=0;
		//
		//  }
		//
		//  for (int i = 0; i < sba_problem->num_observations_; ++i) {
		//    CostFunction* cost_function;
		//    // Each Residual block takes a point and a camera as input and
		//    // outputs a 2 dimensional residual.
		//	int camera_index = sba_problem->camera_index_[i];
		//	double* camera = cameras + camera_block_size * camera_index;   
		//	int point_index = sba_problem->point_index_[i];
		//	int board_id = sba_problem->board_ids()[point_index];
		//	double* board = points + 6 * board_id;
		//
		//    cost_function =
		//    new AutoDiffCostFunction<ReprojectionErrorCalib, 2, 9 ,6>(
		//	new ReprojectionErrorCalib(observations[2 * i + 0], observations[2 * i + 1],sba_problem->init_rot_+sba_problem->camera_index_[i]*4,camera[11],camera[10],sba_problem->board_coords()[point_index*2],sba_problem->board_coords()[point_index*2+1],point_index));
		//		
		//   
		//
		//    // If enabled use Huber's loss function.
		//    LossFunction* loss_function = FLAGS_robustify ? new HuberLoss(1.0) : NULL;
		//
		//    // Each observation correponds to a pair of a camera and a point
		//    // which are identified by camera_index_[i] and point_index_[i]
		//    // respectively.
		//
		//    
		//	
		//	problem->AddResidualBlock(cost_function,
		//    		 	 	 	 	  loss_function,
		//    		 	 	 	 	  camera+1,
		//    		 	 	 	 	  board);
		//	
		//		
		//	
		//	//if(FLAGS_constant_K)problem->SetParameterBlockConstant(camera);
		//	/*if(FLAGS_constant_R)problem->SetParameterBlockConstant(q);*/
		//	if(FLAGS_constant_cameras)problem->SetParameterBlockConstant(camera+1);
		//	if(FLAGS_constant_points)problem->SetParameterBlockConstant(board);	
		//
		//
		//  }
		//  delete [] initrots;
		//
		//}

		void addCalibBoard(BAProblem * sba_problem, Problem* problem){
			double * point_coords = sba_problem->board_coords_;
			int * board_ids = sba_problem->board_ids_;
			double * board_parameters = sba_problem->board_parameters_;
			double* cameras = sba_problem->mutable_cameras();
			double * board_observations = sba_problem->boards_observations_;
			const int camera_block_size = sba_problem->num_camera_parameters_;
			double * initrots = sba_problem->init_rot_;


			for (int i = 0; i < sba_problem->num_boards_observations_; ++i) {
				CostFunction* cost_function;
				// Each Residual block takes a point and a camera as input and
				// outputs a 2 dimensional residual.
				int camera_index = sba_problem->boards_camera_index_[i];
				double* camera = cameras + camera_block_size * camera_index;   
				int point_index = sba_problem->boards_point_index_[i];
				int board_id = sba_problem->board_ids_[point_index];
				double* board = board_parameters + 6 * board_id;

				cost_function =
					new AutoDiffCostFunction<ReprojectionErrorCalib, 2, 9 ,6>(
					new ReprojectionErrorCalib(board_observations[2 * i + 0], board_observations[2 * i + 1],sba_problem->init_rot_+sba_problem->boards_camera_index_[i]*4,camera[11],camera[10],sba_problem->board_coords_[point_index*2],sba_problem->board_coords_[point_index*2+1],point_index));
				// If enabled use Huber's loss function.
				LossFunction* loss_function = FLAGS_robustify ? new HuberLoss(1.0) : NULL;

				// Each observation correponds to a pair of a camera and a point
				// which are identified by camera_index_[i] and point_index_[i]
				// respectively.
				problem->AddResidualBlock(cost_function,
					loss_function,
					camera+1,
					board);
				//if(FLAGS_constant_K)problem->SetParameterBlockConstant(camera);
				/*if(FLAGS_constant_R)problem->SetParameterBlockConstant(q);*/
				if(FLAGS_constant_cameras)problem->SetParameterBlockConstant(camera+1);
				if(FLAGS_constant_points)problem->SetParameterBlockConstant(board);	


			}

		}

		int extractCovarianceBlocks(BAProblem * sba_problem){
			double * pt_ptr;
			double cov_x[3 * 3];
			for (int i = 0; i < sba_problem->num_points_; i++)
			{
				pt_ptr = sba_problem->parameters_ + sba_problem->num_cameras_*sba_problem->num_camera_parameters_+i*3;
				
				if (sba_problem->covariance->GetCovarianceBlock(pt_ptr, pt_ptr, cov_x)){
					sba_problem->cov[sba_problem->num_cameras_*sba_problem->num_camera_parameters_ + i * 3] = cov_x[0];
					sba_problem->cov[sba_problem->num_cameras_*sba_problem->num_camera_parameters_ + i * 3+1] = cov_x[4];
					sba_problem->cov[sba_problem->num_cameras_*sba_problem->num_camera_parameters_ + i * 3+2] = cov_x[8];
				}
			}
			
			double *q, *t, *fx, *ar, *x0y0, *s, *r;
			int offset;
			for (int i = 0; i < sba_problem->num_cameras_; i++)
			{
				double * camera = sba_problem->parameters_ + i*sba_problem->num_camera_parameters_;
				double * res = sba_problem->cov + i*sba_problem->num_camera_parameters_;
				q = camera + 1;
				t = camera + 4;
				fx = camera + 9;
				ar = camera + 12;
				x0y0 = camera + 10;
				s = camera + 13;
				r = camera + 7;

				switch (sba_problem->proj_func_){
				case 0:

					double cov_q[3 * 3];
					double cov_t[3 * 3];
					double cov_fx;
					double cov_ar;
					double cov_x0y0[2 * 2];
					double cov_s;
					double cov_r[2 * 2];

					if (sba_problem->covariance->GetCovarianceBlock(q, q, cov_q)){
						res[1] = cov_q[0];
						res[2] = cov_q[4];
						res[3] = cov_q[8];
					}
					if (sba_problem->covariance->GetCovarianceBlock(t, t, cov_t)){
						res[4] = cov_q[0];
						res[5] = cov_q[4];
						res[6] = cov_q[8];
					}
					if (sba_problem->covariance->GetCovarianceBlock(fx, fx, &cov_fx)){
						res[9] = cov_fx;
					}

					if (sba_problem->covariance->GetCovarianceBlock(ar, ar, &cov_ar)){
						res[12] = cov_ar;
					}

					if (sba_problem->covariance->GetCovarianceBlock(x0y0, x0y0, cov_x0y0)){
						res[10] = cov_x0y0[0];
						res[11] = cov_x0y0[3];
					}
					if (sba_problem->covariance->GetCovarianceBlock(s, s, &cov_s)){
						res[13] = cov_s;
					}
					if (sba_problem->covariance->GetCovarianceBlock(r, r, cov_r)){
						res[7] = cov_r[0];
						res[8] = cov_r[1];
					}

					break;
				case 1:

					offset = 1;
					double cov_res[11 * 11];
					if (sba_problem->covariance->GetCovarianceBlock(camera + offset, camera + offset, cov_res)){

						for (int j = 0; j < 11; j++)
						{
							res[offset + j] = cov_r[j * 12];
						}

					}
					break;
				case 2:
					{
						offset = 1;
						double cov_res[6 * 6];
						if (sba_problem->covariance->GetCovarianceBlock(camera + offset, camera + offset, cov_res)){

							for (int j = 0; j < 6; j++)
							{
								res[offset + j] = cov_r[j * 7];
							}

						}
					}
					break;
				
				case 3:
					{
						offset = 1;
						double cov_res[9 * 9];
						if (sba_problem->covariance->GetCovarianceBlock(camera + offset, camera + offset, cov_res)){

							for (int j = 0; j < 9; j++)
							{
								res[offset + j] = cov_r[j * 10];
							}

						}
					}
					break;
				case 4:
					{
						offset = 1;
						double cov_res[10 * 10];
						if (sba_problem->covariance->GetCovarianceBlock(camera + offset, camera + offset, cov_res)){

							for (int j = 0; j < 10; j++)
							{
								res[offset + j] = cov_r[j * 11];
							}

						}
					}
					break;
				case 5:
					{
						offset = 1;
						double cov_res[7 * 7];
						if (sba_problem->covariance->GetCovarianceBlock(camera + offset, camera + offset, cov_res)){

							for (int j = 0; j < 7; j++)
							{
								res[offset + j] = cov_r[j * 8];
							}

						}
					}
					break;
				case 6:
					{
						offset = 1;
						double cov_res[13 * 13];
						if (sba_problem->covariance->GetCovarianceBlock(camera + offset, camera + offset, cov_res)){

							for (int j = 0; j < 13; j++)
							{
								res[offset + j] = cov_r[j * 14];
							}

						}
					}
					break;
				case 7:
					{
						offset = 1;
						double cov_res[6 * 6];
						if (sba_problem->covariance->GetCovarianceBlock(camera + offset, camera + offset, cov_res)){

							for (int j = 0; j < 6; j++)
							{
								res[offset + j] = cov_r[j * 7];
							}

						}
					}
					break;
				case 8:
					{
						offset = 1;
						double cov_res[11 * 11];
						if (sba_problem->covariance->GetCovarianceBlock(camera + offset, camera + offset, cov_res)){

							for (int j = 0; j < 11; j++)
							{
								res[offset + j] = cov_r[j * 12];
							}

						}
					}
					break;
				case 9:
					{
						offset = 1;
						double cov_res[12 * 12];
						if (sba_problem->covariance->GetCovarianceBlock(camera + offset, camera + offset, cov_res)){

							for (int j = 0; j < 12; j++)
							{
								res[offset + j] = cov_r[j * 13];
							}

						}
					}
					break;
				case 10:
					{
						offset = 1;
						double cov_res[9 * 9];
						if (sba_problem->covariance->GetCovarianceBlock(camera + offset, camera + offset, cov_res)){

							for (int j = 0; j < 9; j++)
							{
								res[offset + j] = cov_r[j * 10];
							}

						}
					}
					break;
				}
				
			}
			
			return 0;
		}
		

		int SolveProblem(BAProblem * ba_problem, Solver::Options options, Solver::Summary * summary){
			Problem problem;
			switch(ba_problem->proj_func_){
			case 0:
				ba_problem->convertT2C();
				ba_problem->convertToSBARadialProblem();
				BuildProblem(ba_problem, &problem);
				break;
			case 1:
			case 2:
			case 3:
			case 4:
			case 5:
				ba_problem->convertToSBAProblem();
				BuildProblem(ba_problem, &problem);
				break;
			case 6:
			case 7:
			case 8:
			case 9:
			case 10:
				ba_problem->convertToSBARadialProblem();
				BuildProblem(ba_problem, &problem);
				break;
			case 12:
				BuildProblem(ba_problem, &problem);
				break;
			case 20:
			case 21:
				ba_problem->convertToBALEaxProblem();
				BuildProblem(ba_problem, &problem);
				break;
			case 22:
			case 23:
			case 24:
			case 25:
				ba_problem->convertToBALQuaternionProblem();
				BuildProblem(ba_problem, &problem);
				break;
			case 30:
			case 31:
				ba_problem->convertToSBAQProblem();
				BuildProblem(ba_problem, &problem);
				break;
			case 40:
			case 41:
			case 42:
				ba_problem->convertToSBARadialProblem();
				BuildProblem(ba_problem, &problem);
				break;
			case 50:
			case 51:
			case 52:
			case 53:
			case 54:
			case 55:
			case 56:
			case 57:
			case 58:
			case 60:
				ba_problem->convertToSBARadialProblem();
				BuildProblem(ba_problem, &problem);
				break;
			case 59:
				ba_problem->convertT2C();
				ba_problem->convertToSBARadialProblem();
				BuildProblem(ba_problem, &problem);
				break;
			case 100:
				BuildProblem(ba_problem, &problem);
				break;
			}


			
			Solve(options, &problem, summary);
			if (ba_problem->covariance_mask!=NULL){
				CHECK(ba_problem->covariance->Compute(ba_problem->covariance_blocks, &problem));
				extractCovarianceBlocks(ba_problem);
			}
			
			switch(ba_problem->proj_func_){
			case 0:
				if (ba_problem->shared_parmask!=NULL)ba_problem->distributeSharedParams();
				ba_problem->convertFromSBARadialProblem();
				ba_problem->convertC2T();
				break;
			case 1:
			case 2:
			case 3:
			case 4:
			case 5:
				ba_problem->convertFromSBAProblem();
				break;
			case 6:
			case 7:
			case 8:
			case 9:
			case 10:
			case 11:
				ba_problem->convertFromSBARadialProblem();
				break;
			case 20:
			case 21:
				ba_problem->convertFromBALEaxProblem();
				break;
			case 22:
			case 23:
			case 24:
			case 25:
				ba_problem->convertFromBALQuaternionProblem();
				break;
			case 30:
			case 31:
				ba_problem->convertFromSBAQProblem();
				break;
			case 40:
			case 41:
			case 42:
				ba_problem->convertFromSBARadialProblem();
				break;
			case 50:
			case 51:
			case 52:
			case 53:
			case 54:
			case 55:
			case 56:
			case 57:
			case 58:
			case 60:
				ba_problem->convertFromSBARadialProblem();
				break;
			case 59: 
				ba_problem->convertC2T();
				ba_problem->convertFromSBARadialProblem();
				break;
			case 100:
				break;
			}

			cout << summary->FullReport() << "\n";
			return 0;


		}

		void SolveProblemFromFile() {
			BAProblem ba_problem;
			ba_problem.proj_func_ = FLAGS_proj_func;
			
			ba_problem.rolling_shutter = FLAGS_rolling_shutter;

			if(FLAGS_estimate_radial){
				ba_problem.loadParametersFromSBARadial(FLAGS_input_cams.c_str(),FLAGS_input_pts.c_str());
			}else{
				ba_problem.loadParametersFromSBANoRadial(FLAGS_input_cams.c_str(),FLAGS_input_pts.c_str());
			}

			if (FLAGS_constant_cameras){
				ba_problem.constant_cameras = new bool[ba_problem.num_cameras_];
				for (int i = 0; i < ba_problem.num_cameras_; i++)
				{
					ba_problem.constant_cameras[i] = true;
				}
			}
			else{
				ba_problem.constant_cameras = new bool[ba_problem.num_cameras_]();
			}
			if (FLAGS_constant_points){
				ba_problem.constant_points = new bool[ba_problem.num_points_];
				for (int i = 0; i < ba_problem.num_points_; i++)
				{
					ba_problem.constant_points[i] = true;
				}
			}
			else{
				ba_problem.constant_points = new bool[ba_problem.num_points_]();
			}
			

			if (FLAGS_rolling_shutter){
				ba_problem.rolling_shutter_orientation = new int[ba_problem.num_cameras_];
				for (int i = 0; i < ba_problem.num_cameras_; i++)
				{
					ba_problem.rolling_shutter_orientation[i] = FLAGS_rs_direction;
				}
			}
			
			/*if(!FLAGS_ref_points.empty()){
				addReferencePoints(&ba_problem,&problem);
			}

			if(!FLAGS_ref_cams.empty()){
				addReferenceCameras(&ba_problem,&problem);
			}*/

			Solver::Options options;
			Solver::Summary summary;
			//check if we are supposed to print out intermediate values
			
			string report;
			SetSolverOptionsFromFlags(&ba_problem, &options);
			options.function_tolerance=1e-60;
			options.parameter_tolerance=1e-60;
			options.gradient_tolerance=1e-60;
			SolveProblem(&ba_problem,options,&summary);			

			if(FLAGS_print_sba){
				if(FLAGS_estimate_radial){
					ba_problem.printSBARadialProblem(FLAGS_output_cams,FLAGS_output_pts);
				}else{
					ba_problem.printSBAProblem(FLAGS_output_cams,FLAGS_output_pts);		
				}
			}  

			if(FLAGS_print_wrl){
				ba_problem.printSBAProblemToWRL(FLAGS_output_wrl);
			}
			if(FLAGS_print_bundler){
				ba_problem.printSBAProblemToBundler(FLAGS_output_bundler);
			}
			if(FLAGS_print_matrices){
				ba_problem.printPMatrices();
			}
			if(FLAGS_print_calib){
				ba_problem.printCalibProblem(FLAGS_output_cams,FLAGS_output_boards);
			}
#ifdef EVAL
			std::ofstream total_time("time.txt");
			total_time << summary.total_time_in_seconds;
			total_time.close();
#endif

		}

		


int SolveProblemFromInterface(BAProblemInterface * problemInterface, std::string * summary_s) {	
	using namespace ceres;

	BAProblem ba_problem;

	if(problemInterface==NULL){
		LOG(ERROR) << "SolveProblem(): failed to process the interface";
		return 1;
	}

	Solver::Options options;
	Solver::Summary summary;
	options.function_tolerance=problemInterface->f_tolerance;
	options.parameter_tolerance=problemInterface->p_tolerance;
	options.gradient_tolerance=problemInterface->g_tolerance;
	options.max_num_iterations = problemInterface->num_iterations;
	options.num_threads = problemInterface->num_threads;
	options.update_state_every_iteration = true;


	ba_problem.loadParametersFromInterface(problemInterface);
	SolveProblem(&ba_problem,options,&summary);
	if (summary_s!=NULL)*summary_s = summary.FullReport();
	
	return 0;



}



void processBundleAdjustmentSettings(BundleAdjustmentSettings * settings){
	//cmp added
	FLAGS_input_cams = settings->input_cams;
	FLAGS_input_pts= settings->input_pts;
	FLAGS_proj_func= settings->proj_func;
	FLAGS_parmask= settings->parmask;
	FLAGS_ref_points= settings->ref_points;
	FLAGS_ref_cams= settings->ref_cams;
	FLAGS_camnames= settings->camnames;
	FLAGS_estimate_radial= settings->estimate_radial;
	FLAGS_print_sba = settings->print_sba;
	FLAGS_print_matrices= settings->print_matrices;
	FLAGS_print_bundler= settings->print_bundler;
	FLAGS_print_wrl= settings->print_wrl;
	FLAGS_output_cams = settings->output_cams;
	FLAGS_output_pts = settings->output_pts;
	FLAGS_output_bundler = settings->output_bundler;
	FLAGS_output_wrl = settings->output_wrl;
	FLAGS_constant_cameras = settings->constant_cameras;
	FLAGS_constant_points = settings->constant_points;

	//ceres original
	FLAGS_trust_region_strategy= settings->trust_region_strategy;
	FLAGS_dogleg= settings->dogleg;

	FLAGS_inner_iterations= settings->inner_iterations;

	FLAGS_blocks_for_inner_iterations= settings->blocks_for_inner_iterations;

	FLAGS_linear_solver= settings->linear_solver;
	FLAGS_preconditioner= settings->preconditioner;
	FLAGS_sparse_linear_algebra_library= settings->sparse_linear_algebra_library;
	FLAGS_ordering= settings->ordering;

	FLAGS_robustify= settings->robustify;
	FLAGS_eta= settings->eta;

	FLAGS_num_threads= settings->num_threads;
	FLAGS_num_iterations= settings->num_iterations;
	FLAGS_max_solver_time= settings->max_solver_time;
	FLAGS_nonmonotonic_steps= settings->nonmonotonic_steps;

	FLAGS_solver_log= settings->solver_log;

}

void initBundleAdjustmentSettings(BundleAdjustmentSettings * settings){
	settings->input_cams = "";
	settings->input_pts = ""; 
	settings->proj_func = 1;
	settings->parmask = "";
	settings->ref_points = ""; 
	settings->ref_cams = ""; 
	settings->camnames = ""; 
	settings->estimate_radial = false;
	settings->print_sba = false;
	settings->print_matrices = false;
	settings->print_bundler = false;
	settings->print_wrl = false;
	settings->output_cams="";
	settings->output_pts="";
	settings->output_bundler="";
	settings->output_wrl="";
	settings->constant_cameras = false;
	settings->constant_points = false;

	//ceres original
	settings->trust_region_strategy = "levenberg_marquardt"; 
	settings->dogleg = "traditional_dogleg";
	settings->inner_iterations = false;
	settings->blocks_for_inner_iterations = "automatic"; 
	settings->linear_solver = "sparse_schur"; 
	settings->preconditioner = "jacobi";
	settings->sparse_linear_algebra_library = "suite_sparse";
	settings->ordering="automatic";
	settings->robustify=false; 

	settings->eta=1e-2; 
	settings->use_block_amd=true; 
	settings->num_threads=1; 
	settings->num_iterations=5; 
	settings->max_solver_time=1e32; 
	settings->nonmonotonic_steps=false;
	settings->solver_log="";
};

BAProblemInterface::BAProblemInterface(){
	parameters = NULL;
	observations = NULL;
	camera_constraints = NULL;
	camera_constraint_weights = NULL;
	error_weights = NULL;
	point_index_ = NULL;
	camera_index_ = NULL;
	n_ref_cams = 0;
	n_ref_pts = 0;
	ref_cam_index = NULL;
	ref_cam_parameters = NULL;
	ref_cam_weights = NULL;
	ref_pts_index = NULL;
	ref_pts_parameters = NULL;
	ref_pts_weights = NULL;
	ncams = 0;
	npts = 0;
	nobs = 0;
	proj_func = 1;
	fixed_parmask = NULL;
	shared_parmask = NULL;
	covariance_mask = NULL;
	num_iterations = 5;
	num_threads = 1;
	p_tolerance = 1e-3;
	g_tolerance = 1e-3;
	f_tolerance = 1e-3;
	robustify = false;
	constant_cameras = NULL;
	constant_points = NULL;
	constrained_cameras = false;
	free_data = false;
	rolling_shutter = false;
	rolling_shutter_orientation = NULL; // 0 for vertical, 1 for horizontal
	im_size = NULL; // image dimensions for calculating the rolling shutter interpolation
	rolling_shutter_frame_time = 0;
	rolling_shutter_capture_time = 0;
	rolling_shutter_v = NULL;
	rolling_shutter_w = NULL;
	rolling_shutter_v_constraints = NULL;
	rolling_shutter_w_constraints = NULL;
	rolling_shutter_v_weight = NULL;
	rolling_shutter_w_weight = NULL;
	q = NULL;
	qRv = NULL;


}

int run_bundle_adjustment() {
	//google::ParseCommandLineFlags(&argc, &argv, true);
	google::InitGoogleLogging("");

	if (FLAGS_input.empty()&&FLAGS_input_cams.empty()) {
		LOG(ERROR) << "Usage: bundle_adjustment_example --input=bal_problem";
		return 1;
	}else{
		CHECK(FLAGS_use_quaternions || !FLAGS_use_local_parameterization)
			<< "--use_local_parameterization can only be used with "
			<< "--use_quaternions.";
		SolveProblemFromFile();
		google::ShutdownGoogleLogging();

	}

	return 0;
}


int main(int argc, char ** argv) {
	google::ParseCommandLineFlags(&argc, &argv, true);
	google::InitGoogleLogging(argv[0]);

	if (FLAGS_input.empty()&&FLAGS_input_cams.empty()) {
		LOG(ERROR) << "Usage: bundle_adjustment_example --input=bal_problem";
		return 1;
	}else{
		CHECK(FLAGS_use_quaternions || !FLAGS_use_local_parameterization)
			<< "--use_local_parameterization can only be used with "
			<< "--use_quaternions.";
		SolveProblemFromFile();

	}

	return 0;
}
