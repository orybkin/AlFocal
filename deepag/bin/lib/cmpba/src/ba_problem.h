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
// Class for loading and holding in memory bundle adjustment problems
// from the BAL (Bundle Adjustment in the Large) dataset from the
// University of Washington.
//
// For more details see http://grail.cs.washington.edu/projects/bal/

#ifndef CERES_EXAMPLES_BAL_PROBLEM_H_
#define CERES_EXAMPLES_BAL_PROBLEM_H_

#include <string>
#include <vector>
#include "cmp_bundle_adjuster.h"
#include <ceres/ceres.h>

class BAProblem {
 public:
 
  explicit BAProblem(){
	  camera_index_=NULL;
	  point_index_=NULL;
	  observations_=NULL;
	  parameters_=NULL;
	  fixed_parmask=NULL;
	  shared_parmask = NULL;
	  covariance_mask = NULL;
	  n_ref_pts_=0;
	  n_ref_cams_=0;
	  proj_func_=1;
	  init_rot_=NULL;
	  nprojs_=NULL;
	  camera_constraints = NULL;
	  error_weights = NULL;
	  num_point_parameters_=3;
	  constant_cameras = NULL;
	  constant_points = NULL;
	  constrained_cameras = false;
	  weighted_errors = false;
	  rolling_shutter = false;
	  rolling_shutter_orientation = NULL;
	  rolling_shutter_frame_time = 0;
	  rolling_shutter_capture_time =0;
	  rolling_shutter_w_constraints=NULL;
	  rolling_shutter_v_constraints = NULL;
	  rolling_shutter_w_weights = NULL;
	  rolling_shutter_v_weights = NULL;
	  rs_fake_last_camera = NULL;
	  im_size=NULL;
	  free_data = 1;
	  cov = NULL;
	  covariance = NULL;
	  cov_options = NULL;
	  q = NULL;
	  qRv = NULL;
	  
  };
  void createBALProblem(const std::string filename, bool use_quaternions);
  void loadParametersFromSBANoRadial(const std::string cameras_filename, const std::string points_filename);
  void loadParametersFromSBARadial(const std::string cameras_filename, const std::string points_filename);
  void loadFixedParmask(const std::string parmask_filename);
  void loadParametersFromInterface(BAProblemInterface * problem);
  void convertToSBAProblem();
  void convertToBALEaxProblem();
  void convertToBALQuaternionProblem();
  void convertToBALQuatProblem();
  void convertToSBARadialProblem();
  void convertFromSBAProblem();
  void convertFromSBARadialProblem();
  void convertFromBALEaxProblem();
  void convertFromBALQuaternionProblem();
  void convertToSBAQProblem();
  void convertToCalibProblem(const std::string boards_filename);
  void convertFromSBAQProblem();
  void convertMulticamToSBAProblem();
  void convertToPancamProblem(std::string parmask_filename);
  void convertToMultiCameraProblem(std::string cameraid_filename);
  void convertPancamToSBAProblem();
  void convertLocalQuatToSBA();
  void convertT2C();
  void convertC2T();
  void distributeSharedParams();
  void normalizeBALProblem();
  void normalizeSBAProblem();
  void printSBAProblem(std::string cam_filename,std::string pts_filename);
  void printSBARadialProblem(std::string cam_filename,std::string pts_filename);
  void printCalibProblem(std::string cam_filename,std::string boards_filename);
  void printSBAProblemToWRL(std::string filename);
  void printSBAProblemToBundler(std::string filename);
  void printPMatrices();
  void printPancamProblem();
  void printRSParameters();
  static void printEvalData(int iteration);
  void loadReferencePoints(std::string filename);
  void loadReferenceCameras(std::string filename);
  void loadBoards(std::string pts_filename, std::string boards_filename);
  ~BAProblem();

  double* mutable_cameras()         { return parameters_;               }
  double* mutable_points()			{    return parameters_  + num_camera_parameters_ * num_cameras_;  }

  int num_cameras_;
  int num_points_;
  int num_observations_;
  int num_parameters_;
  int num_camera_parameters_;
  int num_point_parameters_;
  int num_calibration_parameters_;
  bool use_quaternions_;
  bool robustify;
  bool * constant_cameras;
  bool * constant_points;
  bool constrained_cameras;
  bool weighted_errors;
  bool rolling_shutter;
  int * rolling_shutter_orientation;
  int * im_size;
  double rolling_shutter_frame_time;
  double rolling_shutter_capture_time;
  double * rs_fake_last_camera;
  double * rolling_shutter_v;
  double * rolling_shutter_w;
  double * rolling_shutter_w_constraints;
  double * rolling_shutter_v_constraints;
  double * rolling_shutter_w_weights;
  double * rolling_shutter_v_weights;
  //upvector
  double * qRv;
  double *q;

  bool free_data;

  ceres::Covariance::Options *cov_options;
  ceres::Covariance *covariance;
  std::vector<std::pair<const double*, const double*> > covariance_blocks;
  double * cov;

  double * init_rot_;

  int* point_index_;
  int* camera_index_;
  double* observations_;
  // The parameter vector is laid out as follows
  // [camera_1, ..., camera_n, point_1, ..., point_m]
  double* parameters_;
  int * fixed_parmask;
  int * shared_parmask;
  int * covariance_mask;
  //save numbers of projections for each point to print out SBA file in the end
  int * nprojs_;
  
  double * init_rot_circle;
  double * ref_pts_;
  double * ref_pts_weights_;
  double * ref_cams_;
  double * ref_cams_weights_;

  double * camera_constraints;
  double * camera_constraint_weights;

  double * error_weights;

 
 
  
  int * ref_pts_id_;
  int * ref_cams_id_;
  int n_ref_pts_;
  int n_ref_cams_;
  int proj_func_;

  //boards
  int num_boards_observations_;
  double* boards_observations_;
  int num_boards_;
  int num_board_points_;
  double * board_coords_;
  int * board_ids_;
  double * board_parameters_;
  int num_board_parameters_;

  int* boards_point_index_;
  int* boards_camera_index_;

  //multicam
  int * camids_;
  int * cam_clusters_;
  double * Rs_; //cluster orientations
  double * ts_; //cluster positions
  double * tc_; //shared camera relative positions
  double * Rc_; //shared camera relative orientations
  double * init_rot_cluster_;


  struct {
  	//number of shared parameter blocks
  	int n_blocks;
  	//cluster assigned to each camera in each block
  	int ** block_cam_mask;
  	//pointers to parameters for each camera for each block
  	std::vector<std::vector< double* > > * block_pointers;
  }parmask_;

};



#endif  // CERES_EXAMPLES_BAL_PROBLEM_H_
