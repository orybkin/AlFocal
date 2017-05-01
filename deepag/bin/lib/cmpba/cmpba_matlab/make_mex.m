ceres_bin = ['../../ceres/ceres-bin'];
% ceres_include = ['../../ceres-solver-1.10.0/include']; 

ceres_bin = ['/usr/local/ceres-bin'];
ceres_include = ['/usr/local/ceres-solver-1.12.0/include']; 
cd ../src

command = 'mex GCC="/usr/bin/gcc-4.9" ';
command = 'mex GCC="/usr/bin/gcc-5" ';

files = 'bundle_adjuster_mex.cpp ba_problem.cc bundle_adjuster.cc ';

output = '-output cmp_bundle_adjuster_mex ';

includes = ['-I' ceres_include ' -I/usr/include/gflags -I/usr/include/eigen3 -I' ceres_bin '/config/ '];

lib_dirs = ['-L' ceres_bin '/lib/ '];

libs = '-lceres -lglog -lgflags -lsuitesparseconfig -lcholmod -llapack -lblas -lspqr -lgomp ';


eval([command files output includes lib_dirs libs])

