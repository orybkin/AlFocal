#include <gflags/gflags.h>
#include <glog/logging.h>
#include "ba_problem.h"
#include "snavely_reprojection_error.h"
#include "reprojection_error.h"
#include "quaternion_8d_reprojection_error.h"
#include "PancamReprojectionError.h"
#include "reference3DPointError.h"
#include "referenceCameraError.h"
#include <ceres/ceres.h>
#include <iostream>
#include <fstream>
#include "cmp_bundle_adjuster.h"
#include "camera_constraints_error.h"
#include "RS_reprojection_error.h"
#include "kalibr_radtan_reprojection_error.h"


void BuildProblem(BAProblem* sba_problem, Problem* problem);
int SolveProblem(BAProblem * ba_problem, Solver::Options options, Solver::Summary * summary);