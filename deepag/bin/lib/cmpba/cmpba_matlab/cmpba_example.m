  
load('interface')
addpath('/usr/local/include')
addpath('/usr/local/lib')

%number of cameras
interface.ncams;
%number of points
interface.npts;
%number of observations
interface.nobs;
% PARAMETERS
% N cameras, M points
% [camera1 ... cameraN point1 ... pointM]
% cameraX = [q1 q2 q3 q4 t1 t2 t3 fx x0 y0 ar s r1 r2 r3 r4 r5]
% ar = fy/fx
% r1,...,r5 radial distortion coeeficients, now only r1 and r2 used -
% bundler style
% pointX = [x y z]
interface.parameters;
% next 3 fields for observations follow -
% camera_index,point_index,obsevations
% 3D point (point_index(i)) projects into camera (camera_index(i)) creating
% projection (observations((i-1)*2+1:i*2))
% the order of observation added does not matter, these triplets define
% them
interface.observations;
interface.camera_index;
interface.point_index;
% PROJECTION FUNCTIONS
% defines the function and which parameters to be optimized
% 1-5 standard perspective projection
% 6-10 standard perspective projection with radial
% 1 - K,R,t adjusted
% 2 - K,r fixed, R,t adjsuted
% 3 - ar and s fixed, other adjusted
% 4 - s is fixed, other adjusted
% 5 - ar,s,x0,y0 is fixed, other adjusted
% same holds for 6-10 plus there is radial
% default 1
interface.proj_func = 0;
% for projection function 0 - the 7 camera blocks are [q1,q2,q3][t0,t1,t2][fx][ar][x0,y0][s][r1,r2]
% fixed_parmask - for each block of each camera 0-free/1-fixed
% example of adjusting only camera R and C
fixed_blocks = [0 0 1 1 1 0 0];
for i = 1:interface.ncams
    interface.fixed_parmask((i-1)*7+1:i*7) = fixed_blocks;
end

% which parameter blocks are shared among cameras
% same blocks and order as for fixed parameters
% -1 - block not shared
% 0....N index of a camera with which this parameter block is shared

%first 3 cameras share the internal parameters
shared_blocks = [-1 -1 0 0 0 0 0];
%4th-10th cameras share the internal parameters
shared_blocks2 = [-1 -1 3 3 3 3 3];
for i = 1:3
    interface.shared_parmask((i-1)*7+1:i*7) = shared_blocks;
end
for i = 4:interface.ncams
    interface.shared_parmask((i-1)*7+1:i*7) = shared_blocks2;
end

% array for fixing the cameras
% if constant_cameras(i) = 1 - camera fixed, 0 - camera adjusted, default 0
interface.constant_cameras = zeros(1,interface.ncams);
% array for fixing points, same rules
interface.constant_points = zeros(1,interface.npts);
% maximum number of LM iterations, default 5
interface.num_iterations = 10;
% robustify loss funciton (0 - least squares, 1 - Huber's loss), default 0
interface.robustify = 1;
% number of threads, default 1
interface.num_threads = 4;
interface.verbose = 0;
interface.f_tolerance = 1e-6;
interface.p_tolerance = 1e-6;
interface.g_tolerance = 1e-6;
interface.verbose = 1;
interface.error_weigths = 1*ones(1,interface.nobs);

  
% RUN and obtain adjusted parameters    
[parameters,cost] = cmp_bundle_adjuster_mex(interface);
% this dataset was already adjusted so no update
% lets perturb the points and try again
% interface.parameters(interface.ncams*17+1:end) = interface.parameters(interface.ncams*17+1:end)+ rand(1,interface.npts*3)*0.1;
% [parameters,cost] = cmp_bundle_adjuster_mex(interface);
% parameters(1:17)-interface.parameters(1:17)
