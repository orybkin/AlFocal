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

#include "ba_problem.h"

#include <cstdio>
#include <cstdlib>
#include <string>
#include <glog/logging.h>
#include "ceres/rotation.h"
#include <fstream>
#include <iostream>
#include <Eigen/Eigenvalues>
#include <vector>

#define MAXSTRL 2048





std::string removeExtension(const std::string filename){
    size_t lastdot = filename.find_last_of(".");
    if (lastdot == std::string::npos) return filename;
    return filename.substr(0, lastdot);
}

template<typename T>
void FscanfOrDie(FILE *fptr, const char *format, T *value) {
  int num_scanned = fscanf(fptr, format, value);
  if (num_scanned != 1) {
    LOG(FATAL) << "Invalid UW data file.";
  }
}

template<typename T>
int readNNumbers(FILE *fp, const char *format, T *res, int n){
    for(int i = 0; i <n; i++){
        int j = fscanf(fp,format,res+i);
        if(j==EOF || j != 1) return j;
    }
    return 0;
}

void skipLine(FILE * fp){
    char buf[MAXSTRL];
      while(!feof(fp))
        if(!fgets(buf, MAXSTRL-1, fp) || buf[strlen(buf)-1]=='\n') break;
}

int countLines(std::string filename){
	int cnt=0;
	std::string line;
	std::ifstream in(filename.c_str());
	while(in.good()){
		getline(in, line);
		if(*(line.c_str())!='#'&&line.length())cnt++;
	}
	in.close();
return cnt;
}

void BAProblem::createBALProblem(const std::string filename, bool use_quaternions) {
  FILE* fptr = fopen(filename.c_str(), "r");
  proj_func_=0;
  
  if (!fptr) {
    LOG(FATAL) << "Error: unable to open file " << filename;
    return;
  };

  // This wil die horribly on invalid files. Them's the breaks.
  FscanfOrDie(fptr, "%d", &num_cameras_);
  FscanfOrDie(fptr, "%d", &num_points_);
  FscanfOrDie(fptr, "%d", &num_observations_);

  VLOG(1) << "Header: " << num_cameras_
          << " " << num_points_
          << " " << num_observations_;

  point_index_ = new int[num_observations_];
  camera_index_ = new int[num_observations_];
  observations_ = new double[2 * num_observations_];

  num_parameters_ = 9 * num_cameras_ + 3 * num_points_;
  nprojs_ = new int[num_points_]();
  parameters_ = new double[num_parameters_];

  for (int i = 0; i < num_observations_; ++i) {
    FscanfOrDie(fptr, "%d", camera_index_ + i);
    FscanfOrDie(fptr, "%d", point_index_ + i);
    for (int j = 0; j < 2; ++j) {
      FscanfOrDie(fptr, "%lf", observations_ + 2*i + j);
	}
	nprojs_[point_index_[i]]++;
  }

  for (int i = 0; i < num_parameters_; ++i) {
    FscanfOrDie(fptr, "%lf", parameters_ + i);
  }
  num_camera_parameters_=9;
  use_quaternions_ = use_quaternions;
  if (use_quaternions) {
    // Switch the angle-axis rotations to quaternions.
    num_parameters_ = 10 * num_cameras_ + 3 * num_points_;
    double* quaternion_parameters = new double[num_parameters_];
    double* original_cursor = parameters_;
    double* quaternion_cursor = quaternion_parameters;
    for (int i = 0; i < num_cameras_; ++i) {
		ceres::AngleAxisToQuaternion(original_cursor, quaternion_cursor);
      quaternion_cursor += 4;
      original_cursor += 3;
      for (int j = 4; j < 10; ++j) {
       *quaternion_cursor++ = *original_cursor++;
      }
    }
	num_camera_parameters_=10;
    // Copy the rest of the points.
    for (int i = 0; i < 3 * num_points_; ++i) {
      *quaternion_cursor++ = *original_cursor++;
    }
    // Swap in the quaternion parameters.
    delete []parameters_;
    parameters_ = quaternion_parameters;
  }
  
}

void BAProblem::normalizeBALProblem(){
	for(int i = 0; i<num_observations_;i++){
		int camid = camera_index_[i];
		double * cam = parameters_ + camid*num_camera_parameters_;
		observations_[i*2]/=cam[num_camera_parameters_-3];
		observations_[i*2+1]/=cam[num_camera_parameters_-3];
	}
	for(int i = 0;i<num_cameras_;i++){
		double * cam = parameters_ + i*num_camera_parameters_;
		double f = cam[num_camera_parameters_-3];
		//cam[num_camera_parameters_-2]*=f*f;
		//cam[num_camera_parameters_-1]*=f*f*f*f;
		cam[num_camera_parameters_-3]=1;
	}
}

void BAProblem::normalizeSBAProblem(){
	for(int i = 0; i<num_observations_;i++){
		int camid = camera_index_[i];
		double * cam = parameters_ + camid*num_camera_parameters_;
		observations_[i*2]/=cam[num_camera_parameters_-5];
		observations_[i*2+1]/=cam[num_camera_parameters_-5];
	}
	for(int i = 0;i<num_cameras_;i++){
		double * cam = parameters_ + i*num_camera_parameters_;
		cam[num_camera_parameters_-5]=1;
	}
}

void BAProblem::loadParametersFromSBANoRadial(const std::string cameras_filename, const std::string points_filename) {
  FILE* cfptr = fopen(cameras_filename.c_str(), "r");
  FILE* pfptr = fopen(points_filename.c_str(), "r");
  
  int no_points=0;
  if (!cfptr) {
    LOG(FATAL) << "Error: unable to open camera parameters file " << cameras_filename;
    return;
  };

  if (!pfptr) {
    VLOG(1) << "Did not find file containing points, loading no projections " << points_filename;
	no_points=1;
  };
	
	num_camera_parameters_ = 17;
	int num_camera_file_parameters = 12;

	num_calibration_parameters_=5;
	  


  num_cameras_= countLines(cameras_filename.c_str());
  num_points_ = countLines(points_filename.c_str());
  num_parameters_ = num_camera_parameters_ * num_cameras_ + 3 * num_points_;
  parameters_ = new double[num_parameters_];
  nprojs_ = new int[num_points_];

  if (rolling_shutter){
	  rolling_shutter_v = new double[num_cameras_ * 3];
	  rolling_shutter_w = new double[num_cameras_ * 3];

	  for (int i = 0; i < num_cameras_*3; i++)
	  {
		  rolling_shutter_v[i] = 0;
		  rolling_shutter_w[i] = 0;
	  }
	  
  }


  // Read camera parameters
  int i = 0;
  char ch;
  double temp[12];
  while(!feof(cfptr)&&i<num_cameras_){
    if((ch=fgetc(cfptr))=='#') {
      skipLine(cfptr);
      continue;
    }
    ungetc(ch,cfptr);

    readNNumbers(cfptr,"%lf",temp,num_camera_file_parameters);
	//arrange as q,t,k or q,t,r,k
	//copy t
	memcpy(parameters_ + i*num_camera_parameters_+4,temp+num_calibration_parameters_+4,3*sizeof(double));
	//copy k
	memcpy(parameters_ + i*num_camera_parameters_+7,temp,5*sizeof(double));
	//copy q
    memcpy(parameters_ + i*num_camera_parameters_,temp+num_calibration_parameters_,4*sizeof(double));
	////ar -> fy
	//double fy = parameters_[i*num_camera_parameters_+10]*parameters_[i*num_camera_parameters_+7];
	//parameters_[i*num_camera_parameters_+10]=parameters_[i*num_camera_parameters_+9];
	//parameters_[i*num_camera_parameters_+9]=parameters_[i*num_camera_parameters_+8];
	//parameters_[i*num_camera_parameters_+8]=fy;
	for(int j=num_camera_file_parameters; j<num_camera_parameters_;j++){
		parameters_[i*num_camera_parameters_+j]=0;
	}

    i++;
  }


  double * point_parameters_pointer = parameters_ + num_camera_parameters_*num_cameras_;

  if(no_points){
	  num_points_ = 0;
	  num_observations_=0;  
  }else{
	  // read point parameters and count number of observations
	  num_observations_ = 0;
	  i = 0;
	  int n;
	  while(!feof(pfptr)&&i<num_points_){
		if((ch=fgetc(pfptr))=='#') {
		  skipLine(pfptr);
		  continue;
		}
		ungetc(ch,pfptr);
		readNNumbers(pfptr,"%lf", point_parameters_pointer + i*3,3);
		readNNumbers(pfptr, "%d", &n, 1);
		num_observations_ += n;
		i++;
		skipLine(pfptr);
	  }

	  rewind ( pfptr );

	  observations_ = new double[2 * num_observations_];
	  point_index_ = new int[num_observations_];
	  camera_index_ = new int[num_observations_];

	  i = 0;
	  int curr_obs = 0;
	  double * dummy = new double[3];
	  while(!feof(pfptr)&&i<num_points_){
		if((ch=fgetc(pfptr))=='#') {
		  skipLine(pfptr);
		  continue;
		}
		ungetc(ch,pfptr);
		readNNumbers(pfptr,"%lf", dummy ,3);
		readNNumbers(pfptr, "%d", &n, 1);
		nprojs_[i]=n;
		for(int j=0;j<n;j++){
			point_index_[curr_obs] = i;
			readNNumbers(pfptr,"%d",camera_index_ + curr_obs,1);
			readNNumbers(pfptr,"%lf",observations_ + curr_obs*2,2);
			curr_obs++;
		}
		i++;
		skipLine(pfptr);
	  }

	  delete [] dummy;
  }

  VLOG(1) << "Found: " << num_cameras_
          << " " << num_points_
          << " " << num_observations_;

}

void BAProblem::loadParametersFromSBARadial(const std::string cameras_filename, const std::string points_filename) {
  FILE* cfptr = fopen(cameras_filename.c_str(), "r");
  FILE* pfptr = fopen(points_filename.c_str(), "r");
  int no_points=0;
  if (!cfptr) {
    LOG(FATAL) << "Error: unable to open camera parameters file " << cameras_filename;
    return;
  };

  if (!pfptr) {
    VLOG(1) << "Did not find file containing points, loading no projections " << points_filename;
	no_points=1;
  };
	
	num_camera_parameters_ = 17;

	num_calibration_parameters_=7;
	  


  num_cameras_= countLines(cameras_filename.c_str());
  num_points_ = countLines(points_filename.c_str());
  num_parameters_ = num_camera_parameters_ * num_cameras_ + 3 * num_points_;
  parameters_ = new double[num_parameters_];
  nprojs_ = new int[num_points_];

  // Read camera parameters
  int i = 0;
  char ch;
  double temp[14];
  while(!feof(cfptr)&&i<num_cameras_){
    if((ch=fgetc(cfptr))=='#') {
      skipLine(cfptr);
      continue;
    }
    ungetc(ch,cfptr);

    readNNumbers(cfptr,"%lf",temp,14);
	//arrange as q,t,k or q,t,r,k
	//copy t
	memcpy(parameters_ + i*num_camera_parameters_+4,temp+num_calibration_parameters_+4,3*sizeof(double));
	//copy r
	memcpy(parameters_ + i*num_camera_parameters_+12,temp+5,2*sizeof(double));
	//copy k
	memcpy(parameters_ + i*num_camera_parameters_+7,temp,5*sizeof(double));
	//copy q
    memcpy(parameters_ + i*num_camera_parameters_,temp+num_calibration_parameters_,4*sizeof(double));
	for(int j=14;j<17;j++){
		parameters_[i*num_camera_parameters_+j] = 0;
	}
    i++;
  }


  double * point_parameters_pointer = parameters_ + num_camera_parameters_*num_cameras_;

  if(no_points){
	  num_points_ = 0;
	  num_observations_=0;  
  }else{
	  // read point parameters and count number of observations
	  num_observations_ = 0;
	  i = 0;
	  int n;
	  while(!feof(pfptr)&&i<num_points_){
		if((ch=fgetc(pfptr))=='#') {
		  skipLine(pfptr);
		  continue;
		}
		ungetc(ch,pfptr);
		readNNumbers(pfptr,"%lf", point_parameters_pointer + i*3,3);
		readNNumbers(pfptr, "%d", &n, 1);
		num_observations_ += n;
		i++;
		skipLine(pfptr);
	  }

	  rewind ( pfptr );

	  observations_ = new double[2 * num_observations_];
	  point_index_ = new int[num_observations_];
	  camera_index_ = new int[num_observations_];

	  i = 0;
	  int curr_obs = 0;
	  double * dummy = new double[3];
	  while(!feof(pfptr)&&i<num_points_){
		if((ch=fgetc(pfptr))=='#') {
		  skipLine(pfptr);
		  continue;
		}
		ungetc(ch,pfptr);
		readNNumbers(pfptr,"%lf", dummy ,3);
		readNNumbers(pfptr, "%d", &n, 1);
		nprojs_[i]=n;
		for(int j=0;j<n;j++){
			point_index_[curr_obs] = i;
			readNNumbers(pfptr,"%d",camera_index_ + curr_obs,1);
			readNNumbers(pfptr,"%lf",observations_ + curr_obs*2,2);
			curr_obs++;
		}
		i++;
		skipLine(pfptr);
	  }

	  delete [] dummy;
  }

  VLOG(1) << "Found: " << num_cameras_
          << " " << num_points_
          << " " << num_observations_;

}

void BAProblem::loadFixedParmask(const std::string parmask_filename){
	int i=0;
	fixed_parmask = new int[num_cameras_*7];
	std::ifstream in(parmask_filename.c_str());
	while(i<num_cameras_*7&&!in.eof()){
		in >> fixed_parmask[i]; 
		i++;
	}
	if(i<num_cameras_){
		LOG(ERROR) << "Wrong number of parmask parameters";
	}
}


void BAProblem::distributeSharedParams(){
	for (int i = 0; i < num_cameras_; i++)
	{
		double *camera = parameters_ + i*num_camera_parameters_;
		if (shared_parmask[i * 7]>=0){
			//shared q
			memcpy(camera + 1, parameters_ + shared_parmask[i * 7] * num_camera_parameters_ + 1, 3 * sizeof(double));
		}
		if (shared_parmask[i * 7 + 1] >= 0){
			//shared t
			memcpy(camera + 4, parameters_ + shared_parmask[i * 7+1] * num_camera_parameters_ + 4, 3 * sizeof(double));
		}
		if (shared_parmask[i * 7 + 2] >= 0){
			//shared fx
			memcpy(camera + 9, parameters_ + shared_parmask[i * 7+2] * num_camera_parameters_ + 9,  sizeof(double));
		}
		if (shared_parmask[i * 7 + 3] >= 0){
			//shared ar
			memcpy(camera + 12, parameters_ + shared_parmask[i * 7+3] * num_camera_parameters_ + 12,  sizeof(double));
		}
		if (shared_parmask[i * 7 + 4] >= 0){
			//shared x0y0
			memcpy(camera + 10, parameters_ + shared_parmask[i * 7+4] * num_camera_parameters_ + 10, 2* sizeof(double));
		}
		if (shared_parmask[i * 7 + 5] >= 0){
			//shared s
			memcpy(camera + 13, parameters_ + shared_parmask[i * 7+5] * num_camera_parameters_ + 13, 1 * sizeof(double));
		}
		if (shared_parmask[i * 7 + 6] >= 0){
			//shared r
			memcpy(camera + 7, parameters_ + shared_parmask[i * 7+6] * num_camera_parameters_ + 7, 2 * sizeof(double));
		}
	}
}

void BAProblem::loadParametersFromInterface(BAProblemInterface * problem) {
  proj_func_=problem->proj_func;
  robustify=problem->robustify;
  int no_points=0;
  if (problem->npts==0) {
	no_points=1;
  };

  proj_func_=problem->proj_func;

  if (proj_func_ == 100){
	  num_camera_parameters_ = 21;
  }
  else{
	  num_camera_parameters_ = 17;
  }
	num_calibration_parameters_=10;	  
	 
 
  num_cameras_= problem->ncams;
  num_points_ = problem->npts;
  num_parameters_ = num_camera_parameters_ * num_cameras_ + 3 * num_points_;
  num_observations_=problem->nobs;

  //parameters_ = new double[num_parameters_];
  //observations_ = new double[num_observations_*2];

  //memcpy(parameters_,problem->parameters,num_parameters_*sizeof(double));
  //memcpy(observations_,problem->observations,num_observations_*2*sizeof(double));
  parameters_=problem->parameters;
  observations_=problem->observations;
  camera_index_=problem->camera_index_;
  point_index_=problem->point_index_;
  observations_=problem->observations;
  
  fixed_parmask = problem->fixed_parmask;
  shared_parmask = problem->shared_parmask;
  covariance_mask = problem->covariance_mask;
  if (covariance_mask != NULL){
	  cov = new double[num_cameras_*num_camera_parameters_+num_points_*3];
	  cov_options = new ceres::Covariance::Options();
	  cov_options->algorithm_type = ceres::DENSE_SVD;
	  covariance = new ceres::Covariance(*cov_options);
  }
  if (problem->constant_cameras==NULL){
	  constant_cameras = new bool[problem->ncams](); //default bool = false
  }else{
	  constant_cameras = problem->constant_cameras;
  }
  if (problem->constant_points == NULL){
	  constant_points = new bool[problem->npts](); //default bool = false
  }
  else{
	  constant_points = problem->constant_points;
  }
  free_data = problem->free_data;
  constrained_cameras = problem->constrained_cameras;

  if (problem->error_weights != NULL){
	  weighted_errors = true;
	  error_weights = problem->error_weights;
  }
  else{
	  error_weights = new double[num_observations_];
	  for (int i = 0; i < num_observations_; i++)
	  {
		  error_weights[i] = 1;
	  }
  }

  
  n_ref_cams_ = problem->n_ref_cams;
  n_ref_pts_ = problem->n_ref_pts;
  ref_cams_id_ = problem->ref_cam_index;
  ref_pts_id_ = problem->ref_pts_index;
  ref_cams_ = problem->ref_cam_parameters;
  ref_pts_ = problem->ref_pts_parameters;
  ref_pts_weights_ = problem->ref_pts_weights;
  ref_cams_weights_ = problem->ref_cam_weights;


  rolling_shutter = problem->rolling_shutter;
  rolling_shutter_orientation = problem->rolling_shutter_orientation;
  rolling_shutter_frame_time = problem->rolling_shutter_frame_time;
  rolling_shutter_capture_time = problem->rolling_shutter_capture_time;
  rolling_shutter_v = problem->rolling_shutter_v;
  rolling_shutter_w = problem->rolling_shutter_w;
  rolling_shutter_v_constraints = problem->rolling_shutter_v_constraints;
  rolling_shutter_w_constraints = problem->rolling_shutter_w_constraints;
  rolling_shutter_v_weights = problem->rolling_shutter_v_weight;
  rolling_shutter_w_weights = problem->rolling_shutter_w_weight;

  q = problem->q;
  qRv = problem->q;
  
  im_size = problem->im_size;

  if(constrained_cameras){
	  camera_constraints = problem->camera_constraints;
	  camera_constraint_weights = problem->camera_constraint_weights;  
  }else{
	  camera_constraints=NULL;
	  camera_constraint_weights=NULL;
  }

  double * point_parameters_pointer = parameters_ + num_camera_parameters_*num_cameras_;

  if(no_points){
	  num_points_ = 0;
	  num_observations_=0;  
  }


}

void BAProblem::convertToSBAProblem(){
	init_rot_ = new double[num_cameras_*4];

	for(int i=0;i<num_cameras_;i++){
		//save initial rotation
		memcpy(init_rot_+i*4,parameters_+i*num_camera_parameters_,4*sizeof(double));
		parameters_[i*num_camera_parameters_]=1;
		parameters_[i*num_camera_parameters_+1]=0;
		parameters_[i*num_camera_parameters_+2]=0;
		parameters_[i*num_camera_parameters_+3]=0;	
		//// fx,fy -> ar
		//double ar = parameters_[i*num_camera_parameters_+8]/parameters_[i*num_camera_parameters_+7];
		//parameters_[i*num_camera_parameters_+8]=parameters_[i*num_camera_parameters_+9];
		//parameters_[i*num_camera_parameters_+9]=parameters_[i*num_camera_parameters_+10];
		//parameters_[i*num_camera_parameters_+10]=ar;
	}
}

void BAProblem::convertToSBARadialProblem(){
	init_rot_ = new double[num_cameras_*4];

	for(int i=0;i<num_cameras_;i++){
		//save initial rotation
		memcpy(init_rot_+i*4,parameters_+i*num_camera_parameters_,4*sizeof(double));
		parameters_[i*num_camera_parameters_]=1;
		parameters_[i*num_camera_parameters_+1]=0;
		parameters_[i*num_camera_parameters_+2]=0;
		parameters_[i*num_camera_parameters_+3]=0;	
		// fx,fy -> ar
		/*double ar = parameters_[i*num_camera_parameters_+8]/parameters_[i*num_camera_parameters_+7];
		parameters_[i*num_camera_parameters_+8]=parameters_[i*num_camera_parameters_+9];
		parameters_[i*num_camera_parameters_+9]=parameters_[i*num_camera_parameters_+10];
		parameters_[i*num_camera_parameters_+10]=ar;*/
		//prepare radial
		double r1 = parameters_[i*num_camera_parameters_+12];
		double r2 = parameters_[i*num_camera_parameters_+13];
		double temp[5];
		for(int j=0;j<5;j++){
			temp[j]=parameters_[i*num_camera_parameters_+7+j];
		}
		for(int j=0;j<5;j++){
			parameters_[i*num_camera_parameters_+9+j]=temp[j];
		}

		parameters_[i*num_camera_parameters_ + 7] = r1;
		parameters_[i*num_camera_parameters_ + 8] = r2;
		if (constrained_cameras){
			if (camera_constraints != NULL){
				double cr1 = camera_constraints[i*num_camera_parameters_ + 12];
				double cr2 = camera_constraints[i*num_camera_parameters_ + 13];

				for (int j = 0; j < 5; j++){
					temp[j] = camera_constraints[i*num_camera_parameters_ + 7 + j];
				}
				for (int j = 0; j < 5; j++){
					camera_constraints[i*num_camera_parameters_ + 9 + j] = temp[j];
				}

				camera_constraints[i*num_camera_parameters_ + 7] = cr1;
				camera_constraints[i*num_camera_parameters_ + 8] = cr2;

				double crw1 = camera_constraint_weights[i*num_camera_parameters_ + 12];
				double crw2 = camera_constraint_weights[i*num_camera_parameters_ + 13];

				for (int j = 0; j < 5; j++){
					temp[j] = camera_constraint_weights[i*num_camera_parameters_ + 7 + j];
				}
				for (int j = 0; j < 5; j++){
					camera_constraint_weights[i*num_camera_parameters_ + 9 + j] = temp[j];
				}

				camera_constraint_weights[i*num_camera_parameters_ + 7] = crw1;
				camera_constraint_weights[i*num_camera_parameters_ + 8] = crw2;
			}
		}

	}

	if(rolling_shutter){
		rs_fake_last_camera = new double[num_camera_parameters_];
		memcpy(rs_fake_last_camera,parameters_+(num_cameras_-1)*num_camera_parameters_,num_camera_parameters_*sizeof(double));
	}
}

void BAProblem::convertToBALEaxProblem(){
	double temp[3];
	for (int i = 0; i < num_cameras_; ++i) {
		double * camera = parameters_+num_camera_parameters_*i;
		ceres::QuaternionToAngleAxis(camera, temp);
      for (int j = 0; j < 3; ++j) {
       camera[1+j]=temp[j];
      }
	  //r1,r2 move behind f
	  double r1 = camera[12];
	  double r2 = camera[13];
	  camera[12]=camera[8];
	  camera[13]=camera[9];
	  camera[8]=r1;
	  camera[9]=r2;
    }

}

void BAProblem::convertFromBALEaxProblem(){
	double temp[4];
	for (int i = 0; i < num_cameras_; ++i) {
		double * camera = parameters_+num_camera_parameters_*i;
		ceres::AngleAxisToQuaternion(camera, temp);
      for (int j = 0; j < 4; ++j) {
       camera[j]=temp[j];
      }
	   //r1,r2 move behind s
	  double r1 = camera[8];
	  double r2 = camera[9];
	  camera[8]=camera[12];
	  camera[9]=camera[13];
	  camera[12]=r1;
	  camera[13]=r2;
    }

}

void BAProblem::convertToBALQuaternionProblem(){
	for (int i = 0; i < num_cameras_; ++i) {
		double * camera = parameters_+num_camera_parameters_*i;
	   //r1,r2 move behind s
	  double r1 = camera[12];
	  double r2 = camera[13];
	  camera[12]=camera[8];
	  camera[13]=camera[9];
	  camera[8]=r1;
	  camera[9]=r2;
    }
}

void BAProblem::convertFromBALQuaternionProblem(){
	for (int i = 0; i < num_cameras_; ++i) {
		double * camera = parameters_+num_camera_parameters_*i;
	   //r1,r2 move behind s
	  double r1 = camera[8];
	  double r2 = camera[9];
	  camera[8]=camera[12];
	  camera[9]=camera[13];
	  camera[12]=r1;
	  camera[13]=r2;
    }

}

void BAProblem::convertFromSBAProblem(){
	double tquat[4];
	double localquat[4];
	for(int i =0; i<num_cameras_;i++){
		double * q = parameters_+i*num_camera_parameters_+1;
		localquat[1]= q[0];
		localquat[2]= q[1];
		localquat[3]= q[2];
		localquat[0] = ceres::sqrt(1-q[0]*q[0]-q[1]*q[1]-q[2]*q[2]);
		ceres::QuaternionProduct(localquat, init_rot_ + i * 4, q - 1);
		//ar -> fy
		/*double fy = parameters_[i*num_camera_parameters_+10]*parameters_[i*num_camera_parameters_+7];
		parameters_[i*num_camera_parameters_+10]=parameters_[i*num_camera_parameters_+9];
		parameters_[i*num_camera_parameters_+9]=parameters_[i*num_camera_parameters_+8];
		parameters_[i*num_camera_parameters_+8]=fy;*/
	}	


}

void BAProblem::convertFromSBARadialProblem(){
	double tquat[4];
	double localquat[4];
	for(int i =0; i<num_cameras_;i++){
		double * q = parameters_+i*num_camera_parameters_+1;
		localquat[1]= q[0];
		localquat[2]= q[1];
		localquat[3]= q[2];
		localquat[0] = sqrt(1-q[0]*q[0]-q[1]*q[1]-q[2]*q[2]);
		ceres::QuaternionProduct(localquat, init_rot_ + i * 4, q - 1);
		//prepare radial
		double r1 = parameters_[i*num_camera_parameters_+7];
		double r2 = parameters_[i*num_camera_parameters_+8];
		double temp[5];
		for(int j=0;j<5;j++){
			temp[j]=parameters_[i*num_camera_parameters_+9+j];
		}
		for(int j=0;j<5;j++){
			parameters_[i*num_camera_parameters_+7+j]=temp[j];
		}

		parameters_[i*num_camera_parameters_+12]=r1;
		parameters_[i*num_camera_parameters_+13]=r2;
	}	


}

void BAProblem::convertToSBAQProblem(){
	double * camera,k11;
	double R[9],t[3];
	for(int k =0; k<num_cameras_;k++){
		camera = parameters_+k*num_camera_parameters_;
		k11 = camera[7];
		camera[0]*=sqrt(k11);
		camera[1]*=sqrt(k11);
		camera[2]*=sqrt(k11);
		camera[3]*=sqrt(k11);
		//camera[3]*=k11;
		t[0]=camera[4];
		t[1]=camera[5];
		t[2]=camera[6];
		ceres::QuaternionToRotation(camera, R);
		for(int j=0 ; j<3; ++j){
			  camera[4+j]=-R[j]*t[0]-R[3+j]*t[1]-R[6+j]*t[2];
		}
		//move r1,r2 behind t
		double r1 = camera[12];
		double r2 = camera[13];
		camera[13]=camera[8];
		camera[7]=r1;
		camera[8]=r2;
	}
}



void BAProblem::convertFromSBAQProblem(){
	double * camera,k11;
	double R[9],t[3],*q;
	for(int k =0; k<num_cameras_;k++){
		camera = parameters_+k*num_camera_parameters_;
		q = camera;
		k11 = q[0]*q[0]+q[1]*q[1]+q[2]*q[2]+q[3]*q[3]; 
		q[0]/=sqrt(k11);
		q[1]/=sqrt(k11);
		q[2]/=sqrt(k11);
		q[3]/=sqrt(k11);
		//camera[3]*=k11;
		t[0]=camera[4];
		t[1]=camera[5];
		t[2]=camera[6];
		ceres::QuaternionToRotation(q, R);
		for(int j=0 ; j<3; ++j){
			  camera[4+j]=-R[3*j]*t[0]-R[3*j+1]*t[1]-R[3*j+2]*t[2];
		}
		//move radial back
		double x0=camera[13];
		camera[12]=camera[7];
		camera[13]=camera[8];
		camera[7]=k11;
		camera[8]=x0;
	}
}

void BAProblem::convertToPancamProblem(std::string parmask_filename){
	double R[9];
	double * t;
	int num_pancam_parameters;
	int num_calibration_parameters;
	double * parameters;
	proj_func_=2;
	/*if(estimate_radial_){
		num_pancam_parameters = 20;
		num_calibration_parameters =7;
	}else{
		num_pancam_parameters = 18;
		num_calibration_parameters=5;
	}*/
	parameters = new double[num_cameras_*num_pancam_parameters+num_points_*3];
	//create new parameter vector with suitable size
	for(int k =0; k<num_cameras_;k++){

		double * pancam = parameters+k*num_pancam_parameters;
		double * sba_cam = parameters_+k*num_camera_parameters_;
		register int i,j;
		// copy calibration
		memcpy(pancam,sba_cam+7,num_calibration_parameters);
		//copy orientation
		memcpy(pancam+num_calibration_parameters,sba_cam,4);		  

		// find C = -R'*t;
		ceres::QuaternionToRotation(sba_cam, R);
		  t = sba_cam+4;
		  for(j=0 ; j<3; ++j){
			  pancam[num_calibration_parameters+4+j]=-R[j]*t[0]-R[3+j]*t[1]-R[6+j]*t[2];
		  }
		//set the rest to 0
		  for(i=num_calibration_parameters+7; i<num_pancam_parameters; ++i)
			  pancam[i]=0;
	}
	memcpy(parameters+num_cameras_*num_pancam_parameters,parameters_+num_cameras_*num_camera_parameters_,sizeof(double)*num_points_*3);
	delete [] parameters_;
	parameters_=parameters;
	num_camera_parameters_ = num_pancam_parameters;
	//load the parmask
	std::ifstream parmask(parmask_filename.c_str());
	parmask_.n_blocks =3;
	parmask_.block_cam_mask = new int*[parmask_.n_blocks];
	for (int i = 0; i < parmask_.n_blocks; ++i) {
		parmask_.block_cam_mask[i]= new int[num_cameras_];
	}
	parmask_.block_pointers = new std::vector < std::vector < double * > >[parmask_.n_blocks];
	int cluster=0;

	//over all shared blocks
	for(int i=0;i<parmask_.n_blocks;i++){
		int max_cluster=0;
		//over all cameras
		for(int j=0;j<num_cameras_;j++){
			parmask >> cluster;
			parmask_.block_cam_mask[i][j]=cluster;
			double * cam_parameters = parameters_+j*num_camera_parameters_;
			//found a new cluster
			if(cluster>max_cluster){
				std::vector<double *> temp;
				temp.push_back(cam_parameters);
				parmask_.block_pointers[i].push_back(temp);
				max_cluster=cluster;
			//assign to exisitng cluster
			}else{
				parmask_.block_pointers[i][cluster-1].push_back(cam_parameters);
			}
		}


	}

	//create circles and put cameras on them

		int i,j,k,ncams,dim;
		double *pc, *centers2, rotax[3], rotangle;
		double rotvec[3], S[2][2], Sinv[2][2], C[3], q[4];
		double meanx=0,meany=0,meanz=0,det,r,mean_r,phi;
		double Sxx,Sxy,Syy,Sxxx,Syyy,Sxyy,Syxx,b[2],xy0[2];
		double **centers;
		//process each cluster
		Eigen::Matrix3d cov = Eigen::Matrix3d();

		//first go over all circle centers
		int nclusters = parmask_.block_pointers[0].size();
		for(i =0; i<nclusters; i++){
			meanx=0;
			meany=0;
			meanz=0;
			ncams = parmask_.block_pointers[0][i].size();
			centers = new double*[ncams];
			for(j=0;j<ncams;j++){
				centers[j] = parmask_.block_pointers[0][i][j]+num_calibration_parameters+4;
				meanx+=centers[j][0];
				meany+=centers[j][1];
				meanz+=centers[j][2];
			}
			centers2 = new double[ncams*3];

			meanx=meanx/ncams;
			meany=meany/ncams;
			meanz=meanz/ncams;
			//substract mean
			for(j=0;j<ncams;j++){
				centers[j][0]=centers[j][0]-meanx;
				centers[j][1]=centers[j][1]-meany;
				centers[j][2]=centers[j][2]-meanz;
			}

			//calculate covariance matrix
			for(k=0;k<3;k++){
				cov(k,k) = 0;
				for(j=0;j<ncams;j++){
					cov(k,k)+=centers[j][k]*centers[j][k];
				}
				cov(k,k)=cov(k,k)/ncams;
			}
			cov(0,1)=0;
			for(j=0;j<ncams;j++){
					cov(0,1)+=centers[j][0]*centers[j][1];
			}
			cov(0,1)=cov(1,0)=cov(0,1)/ncams;

			cov(0,2)=0;
			for(j=0;j<ncams;j++){
					cov(0,2)+=centers[j][0]*centers[j][2];
			}
			cov(0,2)=cov(2,0)=cov(0,2)/ncams;

			cov(1,2)=0;
			for(j=0;j<ncams;j++){
					cov(1,2)+=centers[j][1]*centers[j][2];
			}
			cov(1,2)=cov(2,1)=cov(1,2)/ncams;

			dim=3;
			//find eigenvectors and eigenvalues of the covariance matrix
			Eigen::EigenSolver<Eigen::MatrixXd> es;
			es.compute(cov, /* computeEigenvectors = */ true);
			//info = sba_symat_eigvec(covdata, dim, eig);
			Eigen::Vector3d eigenv = (es.eigenvalues()).real();
			double min = eigenv(0);
			int minidx = 0;
			if(min>eigenv(1)){min=eigenv(1); minidx=1;}
			if(min>eigenv(2)){min=eigenv(2); minidx=2;}
			//the norm of the plane corresponds to the eigenvector with the smallest eigenvalue
			Eigen::Matrix3d eigenvec = es.eigenvectors().real();
			Eigen::Vector3d plane_norm = eigenvec.col(minidx);

			//cross product (plane_norm x [0 0 1]) gives us the axis of rotation of the plane
			//cross = {a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]}
			//make sure the norm has z positive
			if(plane_norm(2)<0){
				plane_norm(0)=-plane_norm(0);
				plane_norm(1)=-plane_norm(1);
				plane_norm(2)=-plane_norm(2);
			}
			rotax[0] = plane_norm(1);
			rotax[1] = -plane_norm(0);
			rotax[2] = 0;

			rotangle = acos(plane_norm.data()[2]/sqrt(pow(plane_norm.data()[0],2)+pow(plane_norm.data()[1],2)+pow(plane_norm.data()[2],2)));
			rotvec[0] = rotax[0]*rotangle; rotvec[1] = rotax[1]*rotangle; rotvec[2] = rotax[2]*rotangle;
			//create rotation matrix
			//rotvec2mat(rotvec, &R[0][0]);
			//rotate center's coordiantes
//			for(j=0;j<ncams;j++){
//				centers[j*3] = R[0][0]*centers[j*3]+R[0][1]*centers[j*3+1]+R[0][2]*centers[j*3+2];
//				centers[j*3+1] = R[1][0]*centers[j*3]+R[1][1]*centers[j*3+1]+R[1][2]*centers[j*3+2];
//				centers[j*3+2] = R[2][0]*centers[j*3]+R[2][1]*centers[j*3+1]+R[2][2]*centers[j*3+2];
//			}
			for(j=0;j<ncams;j++){
				ceres::AngleAxisRotatePoint(rotvec, centers[j], centers2 + j * 3);
			}
			Sxx = 0;
			for(j=0;j<ncams;j++){
					Sxx+=centers2[j*3]*centers2[j*3];
			}
			Sxy = 0;
			for(j=0;j<ncams;j++){
					Sxy+=centers2[j*3]*centers2[j*3+1];
			}
			Syy = 0;
			for(j=0;j<ncams;j++){
					Syy+=centers2[j*3+1]*centers2[j*3+1];
			}
			Sxxx = 0;
			for(j=0;j<ncams;j++){
					Sxxx+=centers2[j*3]*centers2[j*3]*centers2[j*3];
			}
			Syyy = 0;
			for(j=0;j<ncams;j++){
					Syyy+=centers2[j*3+1]*centers2[j*3+1]*centers2[j*3+1];
			}
			Sxyy = 0;
			for(j=0;j<ncams;j++){
					Sxyy+=centers2[j*3]*centers2[j*3+1]*centers2[j*3+1];
			}
			Syxx = 0;
			for(j=0;j<ncams;j++){
					Syxx+=centers2[j*3+1]*centers2[j*3]*centers2[j*3];
			}
			S[0][0] = Sxx;
			S[0][1] = Sxy;
			S[1][0] = Sxy;
			S[1][1] = Syy;

			b[0] = (1.0/2)*(Sxxx+Sxyy);
			b[1] = (1.0/2)*(Syyy+Syxx);
			//invert S
			det = S[0][0]*S[1][1]-S[0][1]*S[1][0];
			Sinv[0][0] = S[1][1]/det;
			Sinv[0][1] = -S[0][1]/det;
			Sinv[1][0] = -S[1][0]/det;
			Sinv[1][1] = S[0][0]/det;

			xy0[0]=Sinv[0][0]*b[0]+Sinv[0][1]*b[1];
			xy0[1]=Sinv[1][0]*b[0]+Sinv[1][1]*b[1];

			for(j=0;j<ncams;j++){
				centers2[j*3+0]=centers2[j*3+0]-xy0[0];
				centers2[j*3+1]=centers2[j*3+1]-xy0[1];

			}
			r=sqrt(pow(xy0[0],2) + pow(xy0[1],2)+(Sxx+Syy)/ncams);
			for(j=0;j<3;j++){
				rotvec[j]=-rotvec[j];
			}
			double temp[3];
			//put the circle center back to its orientation
			ceres::AngleAxisRotatePoint(rotvec, xy0, temp);

			C[0]=meanx + temp[0];
			C[1]=meany + temp[1];
			C[2]=meanz + temp[2];
			//save the circle orientation to quaternion
			ceres::AngleAxisToQuaternion(rotvec, q);


			for(j=0;j<ncams;j++){

				pc = parmask_.block_pointers[0][i][j];
				phi = atan2(centers2[j*3+1],centers2[j*3+0]);
				pc[num_calibration_parameters+4] = phi;
				pc[num_calibration_parameters+5] = C[0];
				pc[num_calibration_parameters+6] = C[1];
				pc[num_calibration_parameters+7] = C[2];
				pc[num_calibration_parameters+8] = q[0];
				pc[num_calibration_parameters+9] = q[1];
				pc[num_calibration_parameters+10] = q[2];
				pc[num_calibration_parameters+11] = q[3];
				pc[num_calibration_parameters+12] = r;


			}

			delete [] centers;
			delete [] centers2;

		}
		//average the radius for all clusters that share the same radius

		for(i=0;i<parmask_.block_pointers[2].size();i++){
			mean_r = 0;
			ncams = parmask_.block_pointers[2][i].size();
			for(j=0;j<ncams;j++){
				pc = parmask_.block_pointers[2][i][j];
				mean_r += pc[num_calibration_parameters+12];
			}
			mean_r = mean_r/j;
			for(j=0;j<ncams;j++){

				pc = parmask_.block_pointers[2][i][j];
				pc[num_calibration_parameters+12] = mean_r;
			}
		}

		//save initial orientation of cluster circles
		init_rot_circle = new double[parmask_.block_pointers[1].size()*4];
		for (i = 0; i < parmask_.block_pointers[1].size(); ++i) {
			memcpy(init_rot_circle+i*4,parmask_.block_pointers[1][i][0]+num_calibration_parameters+8,sizeof(double)*4);
		}
}

void BAProblem::convertToMultiCameraProblem(std::string cameraid_filename){
	std::string line;
	std::ifstream in(cameraid_filename.c_str());
	//load camera ids
	camids_ = new int[num_cameras_];
	cam_clusters_ = new int[num_cameras_];
	int ncams_cluster=0;
	for(int i =0;i<num_cameras_;i++){
		in >> camids_[i];
		if(camids_[i]>ncams_cluster)ncams_cluster=camids_[i];
		//meanwhile copy initial rotations back
		memcpy(parameters_+i*num_camera_parameters_,init_rot_+i*4,4*sizeof(double));
	}
	//find full clusters
	int id=1;
	int n_full_clusters = 0;
	int * cluster_beginning= new int[num_cameras_/ncams_cluster];
	int cluster_id=0;
	for(int i=0;i<num_cameras_;i++){
		cam_clusters_[i]=cluster_id;
		if(camids_[i]==id)id++;
		else {
			id=1;
		}
		if(id>ncams_cluster){
		//found a whole cluster!
			cluster_beginning[cluster_id]=i-ncams_cluster+1;
			cluster_id++;
			id=1;
		}
		//correct for offset starting with 1
		camids_[i]--;
	}
	n_full_clusters = cluster_id; 
	//quaternion representations of the cluster orientations
	Rs_ = new double[3*n_full_clusters];
	//cluster centers
	ts_ = new double[3*n_full_clusters];
	// relative camera positions
	tc_ = new double[3*ncams_cluster];
	// relative camera orientations
	Rc_ = new double[3*ncams_cluster];
	//init orientations for clusters
	init_rot_cluster_ = new double[4*n_full_clusters];
	//init orientations for cameras
	init_rot_ = new double[4*ncams_cluster];

	//find cluster representatives - avg rotation and position
	double * cameras = parameters_;
	Eigen::Vector3d avgCs;
	//temporary relative camera positions for calculating the average
	Eigen::Vector3d * avgCc = new Eigen::Vector3d[ncams_cluster];
	Eigen::Vector3d avgEaxs;
	//temporary relative camera positions for calculating the average, eax
	Eigen::Vector3d * avgEaxc = new Eigen::Vector3d[ncams_cluster];
	for(int i=0;i<ncams_cluster;i++){
		avgCc[i] << 0 ,0,0;
		avgEaxc[i]  << 0 ,0,0;
			
	}
	for(int i =0; i<n_full_clusters; i++){
		avgCs << 0,0,0;
		avgEaxs << 0,0,0;
		for(int j=0;j<ncams_cluster;j++){
			double * cam = cameras + (cluster_beginning[i]+j)*num_camera_parameters_;
			double * qc = cam;
			double * t = cam+4;
			double Cc[3];
			double Eaxc[3];
			double Eaxcinv[3];
			//convert quaternion to eax
			ceres::QuaternionToAngleAxis(qc, Eaxc);
			//find Rinv
			Eaxcinv[0]=-Eaxc[0];
			Eaxcinv[1]=-Eaxc[1];
			Eaxcinv[2]=-Eaxc[2];
			//convert t of the camera to C
			ceres::AngleAxisRotatePoint(Eaxcinv, t, Cc);
			
			avgCs(0)+=Cc[0];
			avgCs(1)+=Cc[1];
			avgCs(2)+=Cc[2];
			avgEaxs(0)+=Eaxc[0];
			avgEaxs(1)+=Eaxc[1];
			avgEaxs(2)+=Eaxc[2];

		}
		//center of mass of the cluster
		avgCs/=ncams_cluster;
		avgEaxs/=ncams_cluster;
		
		//compute cluster positions in the cluster orientation frame 
		//and save it to BA parameter array ts
		ceres::AngleAxisRotatePoint(avgEaxs.data(), avgCs.data(), ts_ + 3 * i);
		//save the quaternion representation of cluster orientation to Rs
		ceres::AngleAxisToQuaternion(avgEaxs.data(), init_rot_cluster_ + 4 * i);
		Rs_[3*i]=0;
		Rs_[3*i+1]=0;
		Rs_[3*i+2]=0;
		double Rsinv[4];
		avgEaxs=-avgEaxs;
		ceres::AngleAxisToQuaternion(avgEaxs.data(), Rsinv);


		//for each camera find the relative translation and rotation in the cluster
		for(int j=0;j<ncams_cluster;j++){
			double * cam = cameras + (cluster_beginning[i]+j)*num_camera_parameters_;
			double * qc = cam;
			double * t = cam+4;
			double Cc[3];
			double Eaxc[3];
			double Eaxcinv[3];
			//convert quaternion to eax
			ceres::QuaternionToAngleAxis(qc, Eaxc);
			//find Rinv
			Eaxcinv[0]=-Eaxc[0];
			Eaxcinv[1]=-Eaxc[1];
			Eaxcinv[2]=-Eaxc[2];
			//find Rci = Ri*Rs^T
			ceres::QuaternionProduct(qc, Rsinv, Rc_ + 4 * j);
			//now in eax
			ceres::QuaternionToAngleAxis(Rc_ + 4 * j, Eaxc);
			//convert t of the camera to C
			ceres::AngleAxisRotatePoint(Eaxcinv, t, Cc);
			//from cluster center to camera center in world coord
			Eigen::Vector3d temp;
			temp << Cc[0]-avgCs(0) , Cc[1]-avgCs(1) , Cc[2]-avgCs(2);
			//and now in camera coord
			ceres::AngleAxisRotatePoint(Eaxc, temp.data(), t);
			avgCc[j](0)+=t[0];
			avgCc[j](1)+=t[1];
			avgCc[j](2)+=t[2];
			avgEaxc[j](0)+=Eaxc[0];
			avgEaxc[j](1)+=Eaxc[1];
			avgEaxc[j](2)+=Eaxc[2];
		}
		
		
	}
	//average out the relative camera stuff
	for(int j=0;j<ncams_cluster;j++){
		avgCc[j]/=n_full_clusters;
		avgEaxc[j]/=n_full_clusters;
		ceres::AngleAxisToQuaternion(avgEaxc[j].data(), init_rot_ + 4 * j);
		init_rot_[0];
		init_rot_[1];
		init_rot_[2];
		init_rot_[3];
		Rc_[3*j]=0;
		Rc_[3*j+1]=0;
		Rc_[3*j+2]=0;
		tc_[3*j]=avgCc[j](0);
		tc_[3*j+1]=avgCc[j](1);
		tc_[3*j+2]=avgCc[j](2);
	}



}

void BAProblem::convertToCalibProblem(std::string boardFileName){
	//load boards parameters
	int nboards = countLines(boardFileName);
	std::string line;
	std::ifstream in(boardFileName.c_str());
	int i;
	double * calib_parameters = new double[num_cameras_*num_camera_parameters_+nboards*6];
	int index = num_cameras_*num_camera_parameters_;
	//convert point parameters into board coordinates and board indices
	board_coords_ = new double[num_points_*2];
	board_ids_ = new int[num_points_];
	for(i=0;i<num_points_;i++){
		board_coords_[i*2]=parameters_[index+i*3];
		board_coords_[i*2+1]=parameters_[index+i*3+1];
		board_ids_[i]=(int)(parameters_[index+i*3+2]);
	}

	double angle=0;
	memcpy(calib_parameters,parameters_,sizeof(double)*index);
	for(i=0;i<nboards;i++){
		in >> calib_parameters[index] >> calib_parameters[index+1] >> calib_parameters[index+2];
		in >> angle;
		calib_parameters[index]*=angle;
		calib_parameters[index+1]*=angle;
		calib_parameters[index+2]*=angle;		
		in >> calib_parameters[index+3];
		in >> calib_parameters[index+4];
		in >> calib_parameters[index+5];
		index+=6;
	}
	delete [] parameters_;
	parameters_ = calib_parameters;
	num_boards_=nboards;
}

void BAProblem::convertPancamToSBAProblem(){
	double * sba_parameters;
	int num_sba_cam_parameters;
	int num_calib_parameters;

	/*if(estimate_radial_){
		num_sba_cam_parameters = 14;
		num_calib_parameters = 7;
	}else{
		num_sba_cam_parameters = 12;
		num_calib_parameters = 5;
	}*/

	sba_parameters = new double[num_cameras_*num_sba_cam_parameters+num_points_*3];
	memcpy(sba_parameters+num_cameras_*num_sba_cam_parameters,parameters_+num_cameras_*num_camera_parameters_,sizeof(double)*num_points_*3);
	for(int i=0;i<num_cameras_;i++){
		double * sba_cam = sba_parameters+i*num_sba_cam_parameters;
		double * pancam = parameters_+i*num_camera_parameters_;
		//copy calibration
		for (int j = 0; j < num_calib_parameters; ++j){
			sba_cam[j] = pancam[j];
		}

		//compute orientation of the camera
		double tquat[4];
		double localquat[4];
		double * q =&(parameters_[i*num_camera_parameters_+num_calib_parameters+1]);
		localquat[1]= q[0];
		localquat[2]= q[1];
		localquat[3]= q[2];
		localquat[0] = sqrt(1-q[0]*q[0]-q[1]*q[1]-q[2]*q[2]);
		ceres::QuaternionProduct(init_rot_ + i * 4, localquat, tquat);

		//save the orientation for later
		const double quat[4] = {tquat[0],tquat[1],tquat[2],tquat[3]};

		//copy only the local orientation, the total will be computed in printSBAProblem()
		for (int j = num_calib_parameters; j < num_calib_parameters+4; ++j) {
			sba_cam[j]=pancam[j];
		}

		//compute orientation of the circle
		int cluster_index = parmask_.block_cam_mask[1][i]-1;
		q =parmask_.block_pointers[1][cluster_index][0]+num_calib_parameters+9;
		localquat[1]= q[0];
		localquat[2]= q[1];
		localquat[3]= q[2];
		localquat[0] = sqrt(1-q[0]*q[0]-q[1]*q[1]-q[2]*q[2]);
		ceres::QuaternionProduct(init_rot_circle + (cluster_index)* 4, localquat, tquat);

		const double quatc[4] = {tquat[0],tquat[1],tquat[2],tquat[3]};

		//calculate C of the camera
		double phi = pancam[num_calib_parameters+4];
		cluster_index = parmask_.block_cam_mask[0][i]-1;
		double * Cc = parmask_.block_pointers[0][cluster_index][0]+num_calib_parameters+5;
		cluster_index = parmask_.block_cam_mask[2][i]-1;
		double r = parmask_.block_pointers[2][cluster_index][0][num_calib_parameters+12];
		double Ctemp[3];
		const double temp[3] = {cos(phi)*r,sin(phi)*r,0};


		ceres::QuaternionRotatePoint(quatc, temp, Ctemp);
		const double C[]={Ctemp[0]+Cc[0],Ctemp[1]+Cc[1],Ctemp[2]+Cc[2]};
		double t[3];
		ceres::QuaternionRotatePoint(quat, C, t);
		t[0]=-t[0];
		t[1]=-t[1];
		t[2]=-t[2];
		for (int j = 0; j < 3; ++j) {
			sba_cam[num_sba_cam_parameters-3+j] = t[j];
		}

	}
	free(parameters_);
	parameters_ = sba_parameters;
	num_camera_parameters_ = num_sba_cam_parameters;

}
void BAProblem::convertMulticamToSBAProblem(){
	//calculate position and orientation for each camera

	int clusterid,cam_clusterid;
	double * camera,*qc,*tc,*qs,*ts,*qc0,*qs0;

	for(int i=0;i<num_cameras_;i++){
		clusterid =cam_clusters_[i];
		cam_clusterid = camids_[i];
		camera = parameters_ + num_camera_parameters_ * i;
		qc=Rc_+cam_clusterid*3;
		tc=tc_+cam_clusterid*3;
		qs=Rs_+clusterid*3;
		ts=ts_+clusterid*3;
		qc0 = init_rot_+cam_clusterid*4;
		qs0 = init_rot_cluster_+clusterid*4;
		double w = sqrt(1-qc[0]*qc[0]-qc[1]*qc[1]-qc[2]*qc[2]);
		double qclocal[4]= {w,qc[0],qc[1],qc[2]};
		w = sqrt(1-qs[0]*qs[0]-qs[1]*qs[1]-qs[2]*qs[2]);
		double qslocal[4]= {w,qs[0],qs[1],qs[2]};
		double qctotal[4];
		ceres::QuaternionProduct(qc0, qclocal, qctotal);
		double qstotal[4];
		ceres::QuaternionProduct(qs0, qslocal, qstotal);
		double qtotal[4];
		ceres::QuaternionProduct(qctotal, qstotal, qtotal);
		double temp[3];
		ceres::QuaternionRotatePoint(qctotal, ts, temp);
		double t[3];
		t[0]=temp[0]+tc[0];
		t[1]=temp[1]+tc[1];
		t[2]=temp[2]+tc[2];
		memcpy(camera,qtotal,4*sizeof(double));
		memcpy(camera+4,t,3*sizeof(double));

	}

}

void BAProblem::convertT2C(){
	double * camera;
	double R[9], t[3];
	for (int k = 0; k < num_cameras_; k++){
		camera = parameters_ + num_camera_parameters_ * k;
		t[0] = camera[4];
		t[1] = camera[5];
		t[2] = camera[6];
		ceres::QuaternionToRotation(camera, R);
		for (int j = 0; j < 3; ++j){
			camera[4 + j] = -R[j] * t[0] - R[3 + j] * t[1] - R[6 + j] * t[2];
		}
	}

}

void BAProblem::convertC2T(){
	double * camera;
	double R[9], C[3];
	for (int k = 0; k < num_cameras_; k++){
		camera = parameters_ + num_camera_parameters_ * k;
		C[0] = camera[4];
		C[1] = camera[5];
		C[2] = camera[6];
		ceres::QuaternionToRotation(camera, R);
		for (int j = 0; j < 3; ++j){
			camera[4 + j] = -R[j*3] * C[0] - R[j*3+1] * C[1] - R[j*3+2] * C[2];
		}
	}

}


void BAProblem::printSBAProblem(std::string cam_filename,std::string pts_filename){
	std::ofstream outCams(cam_filename.c_str());
	std::ofstream outPts(pts_filename.c_str());
	outCams.setf(std::ios::scientific);
	outPts.setf(std::ios::scientific);
	outCams.precision(8);
	outPts.precision(8);
	int num_calibration_parameters;
	/*if(estimate_radial_){
		num_calibration_parameters=7;
	}else{
		num_calibration_parameters=5;
	}*/

	//print cameras
	for(int i =0; i<num_cameras_;i++){
		double * camera = parameters_+i*num_camera_parameters_;
		for(int j=0;j<5;j++){
			outCams << camera[7+j] << " ";
		}
		for(int j=0;j<6;j++){
			outCams << camera[j]  << " ";
		}

		outCams << camera[6] << "\n";
	}
	outCams.close();
	//print points
	double * pointspar = parameters_ + num_cameras_*num_camera_parameters_;
	int obs =0;
	for(int i =0; i<num_points_;i++){
		//x y z
		outPts << pointspar[i*3] << " " << pointspar[i*3+1] << " " << pointspar[i*3+2] << " ";
		//number of projections
		outPts << nprojs_[i] << " ";
		// projections
		for(int j =0; j<nprojs_[i];j++){
			//camera id
			outPts << camera_index_[obs] << " ";
			//u,v
			outPts << observations_[obs*2] << " " << observations_[obs*2+1] << " ";
			obs++;
		}
		outPts << "\n";
	}
	outPts.close();

}

void BAProblem::printSBARadialProblem(std::string cam_filename,std::string pts_filename){
	std::ofstream outCams(cam_filename.c_str());
	std::ofstream outPts(pts_filename.c_str());
	outCams.setf(std::ios::scientific);
	outPts.setf(std::ios::scientific);
	outCams.precision(8);
	outPts.precision(8);
	int num_calibration_parameters;
	/*if(estimate_radial_){
		num_calibration_parameters=7;
	}else{
		num_calibration_parameters=5;
	}*/

	//print cameras
	for(int i =0; i<num_cameras_;i++){
		double * camera = parameters_+i*num_camera_parameters_;
		for(int j=0;j<7;j++){
			outCams << camera[7+j] << " ";
		}
		for(int j=0;j<6;j++){
			outCams << camera[j]  << " ";
		}

		outCams << camera[6] << "\n";
	}
	outCams.close();
	//print points
	double * pointspar = parameters_ + num_cameras_*num_camera_parameters_;
	int obs =0;
	for(int i =0; i<num_points_;i++){
		//x y z
		outPts << pointspar[i*3] << " " << pointspar[i*3+1] << " " << pointspar[i*3+2] << " ";
		//number of projections
		outPts << nprojs_[i] << " ";
		// projections
		for(int j =0; j<nprojs_[i];j++){
			//camera id
			outPts << camera_index_[obs] << " ";
			//u,v
			outPts << observations_[obs*2] << " " << observations_[obs*2+1] << " ";
			obs++;
		}
		outPts << "\n";
	}
	outPts.close();

}

void BAProblem::printCalibProblem(std::string cam_filename,std::string boards_filename){
	std::ofstream outCams(cam_filename.c_str());
	std::ofstream outPts(boards_filename.c_str());
	outCams.setf(std::ios::scientific);
	outPts.setf(std::ios::scientific);
	outCams.precision(8);
	outPts.precision(8);
	int num_calibration_parameters;
	/*if(estimate_radial_){
		num_calibration_parameters=7;
	}else{
		num_calibration_parameters=5;
	}*/

	//print cameras
	for(int i =0; i<num_cameras_;i++){
		double tquat[4];
		double localquat[4];
		double temp[7];
		//save the q,t
		memcpy(temp,parameters_+i*num_camera_parameters_,7*sizeof(double));
		//copy k to the front
		memcpy(parameters_+i*num_camera_parameters_, parameters_+i*num_camera_parameters_+7,num_calibration_parameters_*sizeof(double));
		//put q,t in the back
		memcpy(parameters_+i*num_camera_parameters_+num_calibration_parameters_,temp,7*sizeof(double));
		double * q =&(parameters_[i*num_camera_parameters_+num_calibration_parameters+1]);
		localquat[1]= q[0];
		localquat[2]= q[1];
		localquat[3]= q[2];
		localquat[0] = sqrt(1-q[0]*q[0]-q[1]*q[1]-q[2]*q[2]);
		ceres::QuaternionProduct(init_rot_ + i * 4, localquat, tquat);
		q = &(parameters_[i*num_camera_parameters_+num_calibration_parameters]);
		for(int j=0;j<num_camera_parameters_-1;j++){
			if(j>(num_calibration_parameters-1)&&j<(num_calibration_parameters+4)){
				outCams << tquat[j-num_calibration_parameters] << " ";
				q[j-num_calibration_parameters] = tquat[j-num_calibration_parameters];
			}else{
				outCams << parameters_[i*num_camera_parameters_+j] << " ";
			}
		}
		outCams << parameters_[i*num_camera_parameters_+num_camera_parameters_-1] << "\n";

	}
	outCams.close();
	//print boards
	double * boardspar = board_parameters_;
	double eax[4];
	for(int i =0; i<num_boards_;i++){
		//eax
		double norm = sqrt(boardspar[i*6]*boardspar[i*6]+boardspar[i*6+1]*boardspar[i*6+1]+boardspar[i*6+2]*boardspar[i*6+2]); 
		outPts << boardspar[i*6]/norm << " " << boardspar[i*6+1]/norm << " " << boardspar[i*6+2]/norm << " " << norm << " ";
		//t
		outPts << boardspar[i*6+3] << " " << boardspar[i*6+4] << " " << boardspar[i*6+5] << " ";
		outPts << "\n";
	}
	outPts.close();
}

void BAProblem::printSBAProblemToBundler(std::string filename){
	std::ofstream out(filename.c_str());
	int num_calib_parameters;

	Eigen::Matrix<double,3,3,Eigen::RowMajor> R;
	Eigen::Vector3d t;
	Eigen::Matrix<double,3,3,Eigen::RowMajor> rotateX;
	rotateX << 1, 0, 0, 0, -1, 0, 0, 0, -1;
	out.setf(std::ios::scientific,std::ios::floatfield);
	out.precision(8);
	out << "#BUndler file v0.3\n";
	out << num_cameras_ << " " << num_points_ << "\n";

	for(int i=0;i<num_cameras_;i++){
		double * camera = parameters_+i*num_camera_parameters_;
		out << camera[0] << " ";
		/*if(estimate_radial_){
			out << camera[5] << " " << camera[6];
			num_calib_parameters = 7;
		}else{
			out << 0 << " " << 0;
			num_calib_parameters =5;
		}*/
		out << "\n";
		double * q = camera+num_calib_parameters;
		double * T = camera+num_calib_parameters+4;
		Eigen::Vector3d t(T[0],T[1],T[2]);

		ceres::QuaternionToRotation(q, R.data());
		R = rotateX*R;
		t = rotateX*t;

		out << R(0,0) << " " << R(0,1) << " " << R(0,2) << "\n";
		out << R(1,0) << " " << R(1,1) << " " << R(1,2) << "\n";
		out << R(2,0) << " " << R(2,1) << " " << R(2,2) << "\n";
		out << t(0) << " " << t(1) << " " << t(2) << "\n";

	}

	double * pointspar = parameters_ + num_cameras_*num_camera_parameters_;
	int obs =0;
		for(int i =0; i<num_points_;i++){
			//x y z
			out << pointspar[i*3] << " " << pointspar[i*3+1] << " " << pointspar[i*3+2] << "\n";
			//color
			out << 255/2 << " " << 255/2 << " " << 255/2 << "\n";
			//number of projections
			out << nprojs_[i] << " ";
			// projections
			for(int j =0; j<nprojs_[i];j++){
				//camera id
				out << camera_index_[obs] << " ";
				//u,v
				out << observations_[obs*2] << " " << observations_[obs*2+1] << " ";
				obs++;
			}
			out << "\n";
		}

	out.close();
}

void BAProblem::printPMatrices(){
	int num_calib_parameters;

	Eigen::Matrix<double,3,3,Eigen::RowMajor> R;
	Eigen::Matrix<double,3,3,Eigen::RowMajor> K;
	Eigen::Matrix<double,3,4,Eigen::RowMajor> P;
	Eigen::Vector3d t;
	//Eigen::Matrix<double,3,3,Eigen::RowMajor> rotateX;
	//rotateX << 1, 0, 0, 0, -1, 0, 0, 0, -1;

//	std::ifstream in(camnames.c_str());

	for(int i=0;i<num_cameras_;i++){
			char p_name[50];
//			std::string pic_name;
//			in >> pic_name;
//			pic_name = removeExtension(pic_name);
			sprintf(p_name,"%05d_P.txt",i+1);
			double * camera = parameters_+i*num_camera_parameters_;
			if(camera[0]!=0){
				std::ofstream out(p_name);
				out.setf(std::ios::scientific,std::ios::floatfield);
				out.precision(8);
				/*if(estimate_radial_){
					num_calib_parameters = 7;
				}else{
					num_calib_parameters =5;
				}*/
				double * q = camera+num_calib_parameters;
				double * T = camera+num_calib_parameters+4;
				Eigen::Vector3d t(T[0],T[1],T[2]);

				ceres::QuaternionToRotation(q, R.data());

				//R = rotateX*R;
				//t = rotateX*t;
				K << camera[0],0,camera[1],
				     0, camera[3]*camera[0], camera[2],
				     0,0,1;


				P.topLeftCorner(3,3) = K*R;
				P.bottomRightCorner(3,1) = K*t;

				//DEBUG OF P MATRICES
//
//				std::cout << "K: \n";
//
//				std::cout << K(0,0) << " " << K(0,1) << " " << K(0,2)  << " \n";
//				std::cout << K(1,0) << " " << K(1,1) << " " << K(1,2)  << " \n";
//				std::cout << K(2,0) << " " << K(2,1) << " " << K(2,2)  << " \n";
//
//				std::cout << "R: \n";
//
//				std::cout << R(0,0) << " " << R(0,1) << " " << R(0,2)  << " \n";
//				std::cout << R(1,0) << " " << R(1,1) << " " << R(1,2)  << " \n";
//				std::cout << R(2,0) << " " << R(2,1) << " " << R(2,2)  << " \n";
//
//				Eigen::Vector3d C;
//				C = -R.transpose()*t;
//
//				std::cout << "C: \n";
//
//				std::cout << C(0) << " " << C(1) << " " << C(2)  << " \n";


				out << P(0,0) << " " << P(0,1) << " " << P(0,2) << " " << P(0,3) << " \n";
				out << P(1,0) << " " << P(1,1) << " " << P(1,2) << " " << P(1,3) << " \n";
				out << P(2,0) << " " << P(2,1) << " " << P(2,2) << " " << P(2,3) << " \n";
				out.close();
			}

		}
}

void BAProblem::printRSParameters(){

}

void BAProblem::printSBAProblemToWRL(std::string filename)
{
	FILE *f = fopen(filename.c_str(),"w");
	fprintf(f, "#VRML V2.0 utf8\n");
	fprintf(f, "Group {\nchildren [\n");
	fprintf(f, "\tWorldInfo {\n");
	fprintf(f, "\tinfo \"test\"\n");
	fprintf(f, "\t}\n\n");
	fprintf(f, "Background {\nskyColor 1 1 1 \n}\n");


	fprintf(f, "Background {\nskyColor 1 1 1 \n}\n");
	fprintf(f, "\tTransform {\n");
	fprintf(f, "\t\ttranslation %g %g %g\n", 0.0, 0.0, 0.0);
	fprintf(f, "\t\tscale %g %g %g\n", 1.0, 1.0, 1.0);
	fprintf(f, "\n\t\tchildren [ Shape {\n");
	fprintf(f, "\n\t\tappearance Appearance { material Material {emissiveColor 1 1 1 }}\n");
	fprintf(f, "\t\t\tgeometry PointSet {\n");

	fprintf(f, "\t\t\t\tcoord Coordinate { point [\n");


	//printf("writing pts\n");
	double * points = mutable_points();
	for (int i=0;i<num_points_;i++)
	{
		fprintf(f, "\t\t\t\t\t%7.8f %7.8f %7.8f \n", points[i*3], points[i*3+1], points[i*3+2]);

	};
	fprintf(f, "\t\t\t\t]}\n\n");
	fprintf(f, "\t\t\t\tcolor Color { color [\n");
	//printf("writing cls\n");
	//load all pts for reference camera
	for (int i=0;i<num_points_;i++)
		{
			fprintf(f, "\t\t\t\t\t%7.8f %7.8f %7.8f \n", points[i*3], points[i*3+1], points[i*3+2]);

		};

	fprintf(f, "\t\t\t\t]}\n\n");

	fprintf(f, "\t\t\t} # end of PointSet\n\n");

	fprintf(f,"\t\t}]#end of children shape\n\t}#end of transform\n] # end of children\n} # end of group\n");


	//printfGroupCameras(f,mp, meanpixsize / 10.0, shift, step);
	//printfGroupCameras(f,mp, 0.01, 0, 1);



	//printf("npts: %i\n",pts->size());
	fclose(f);

}


void BAProblem::loadReferencePoints(std::string filename){
	int npts;
	npts = countLines(filename);
	std::ifstream in(filename.c_str());
	ref_pts_= new double[npts*3];
	ref_pts_weights_ = new double[npts];
	ref_pts_id_ = new int[npts];
	for (int i = 0; i < npts; ++i) {
		in >> ref_pts_id_[i] >> ref_pts_weights_[i];
		in >> ref_pts_[i*3] >> ref_pts_[i*3+1] >> ref_pts_[i*3+2];
	}
	n_ref_pts_ = npts;
}

void BAProblem::loadReferenceCameras(std::string filename){
	int ncams;
	ncams= countLines(filename);
	std::ifstream in(filename.c_str());
	ref_cams_= new double[ncams*3];
	ref_cams_weights_ = new double[ncams];
	ref_cams_id_ = new int[ncams];
	double * q;
	double temp[3];
	int num_calibration_parameters;
		/*if(estimate_radial_){
			num_calibration_parameters=7;
		}else{
			num_calibration_parameters=5;
		}*/
	for (int i = 0; i < ncams; ++i) {
		in >> ref_cams_id_[i] >> ref_cams_weights_[i];
		//q = parameters_+num_camera_parameters_*ref_cams_id_[i]+num_calibration_parameters;
		//const double cam_orientation[4] = {q[0],q[1],q[2],q[3]};
		in >> ref_cams_[i*3] >> ref_cams_[i*3+1] >> ref_cams_[i*3+2];
		//const double C[3] = {-temp[0], -temp[1] ,-temp[2]};
		//QuaternionRotatePoint(cam_orientation,C,&(ref_cams_[i*3]));
	}
	n_ref_cams_ = ncams;
}

void BAProblem::loadBoards(std::string pts_filename,std::string boards_filename){
	num_board_points_ = countLines(pts_filename.c_str());
	num_boards_ = countLines(boards_filename.c_str());
	std::string line;
	std::ifstream in(boards_filename.c_str(), std::ifstream::in);
	int i;
	int index=0;
	double angle=0;
	board_parameters_ = new double[num_boards_*6];
	for(i=0;i<num_boards_;i++){
		in >> board_parameters_[index] >> board_parameters_[index+1] >> board_parameters_[index+2];
		in >> angle;
		board_parameters_[index]*=angle;
		board_parameters_[index+1]*=angle;
		board_parameters_[index+2]*=angle;		
		in >> board_parameters_[index+3];
		in >> board_parameters_[index+4];
		in >> board_parameters_[index+5];
		index+=6;
	}
	in.close();


  // read point parameters and count number of observations
  FILE* pfptr = fopen(pts_filename.c_str(), "r");
  num_boards_observations_ = 0;
  i = 0;
  int n;
  char ch;
  board_coords_ = new double[num_board_points_*2];
	board_ids_ = new int[num_board_points_];
  while(!feof(pfptr)&&i<num_board_points_){
    if((ch=fgetc(pfptr))=='#') {
      skipLine(pfptr);
      continue;
    }
    ungetc(ch,pfptr);
	readNNumbers(pfptr,"%lf", board_coords_ + i*2,2);
	readNNumbers(pfptr,"%d", board_ids_+i,1);
    readNNumbers(pfptr, "%d", &n, 1);
    num_boards_observations_ += n;
    i++;
    skipLine(pfptr);
  }

  rewind ( pfptr );

  boards_observations_ = new double[2 * num_boards_observations_];
  boards_point_index_ = new int[num_boards_observations_];
  boards_camera_index_ = new int[num_boards_observations_];

  i = 0;
  int curr_obs = 0;
  double * dummy = new double[3];
  while(!feof(pfptr)&&i<num_board_points_){
    if((ch=fgetc(pfptr))=='#') {
      skipLine(pfptr);
      continue;
    }
    ungetc(ch,pfptr);
    readNNumbers(pfptr,"%lf", dummy ,3);
    readNNumbers(pfptr, "%d", &n, 1);
    //nprojs_[i]=n;
    for(int j=0;j<n;j++){
        boards_point_index_[curr_obs] = i;
        readNNumbers(pfptr,"%d",boards_camera_index_ + curr_obs,1);
		readNNumbers(pfptr,"%lf",boards_observations_ + curr_obs*2,2);
        curr_obs++;
    }
    i++;
    skipLine(pfptr);
  }
  
  delete [] dummy;
}

BAProblem::~BAProblem() {
  if(free_data){
	  if(point_index_!=NULL)delete []point_index_;
	  if(camera_index_!=NULL)delete []camera_index_;
	  if(observations_!=NULL)delete []observations_;
	  if(parameters_!=NULL)delete []parameters_;
	  if(nprojs_!=NULL)delete []nprojs_;
	  if(fixed_parmask!=NULL)delete []fixed_parmask;
	  delete[] constant_cameras;
	  delete[] constant_points;
  }
  
  if(init_rot_!=NULL){
	delete [] init_rot_;
  }
  if(n_ref_pts_>0){
	  delete [] ref_pts_;
	  delete [] ref_pts_id_;
	  delete [] ref_pts_weights_;
  }

  if(n_ref_cams_>0){
	  delete [] ref_cams_;
	  delete [] ref_cams_id_;
	  delete [] ref_cams_weights_;
  }

}

