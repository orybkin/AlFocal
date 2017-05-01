ceres_bin = ['path_to_ceres_bin_folder'];
ceres_include = ['path_to_ceres_solver/include']; 

command = 'mex ';

files = 'bundle_adjuster_mex.cpp ba_problem.cc bundle_adjuster.cc ';

output = '-output cmp_bundle_adjuster_mex ';

includes = ['-I' ceres_include ' -I/usr/include/gflags -I/usr/include/eigen3 -I' ceres_bin '/config/ '];

lib_dirs = ['-L' ceres_bin '/lib/ '];

libs = '-lceres -lglog -lgflags -lsuitesparseconfig -lcholmod -llapack -lblas -lspqr -lgomp ';

eval([command files output includes lib_dirs libs])

