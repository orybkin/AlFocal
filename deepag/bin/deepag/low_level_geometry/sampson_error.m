function [ err ] = sampson_error( F,u1,u2,mth )
%SAMPSON_ERROR Summary of this function goes here
%   Detailed explanation goes here
err=sampson(F,u1,u2);
if nargin<4
    mth='^2avg';
end
switch mth
    case '^2avg'
        err=norm(err)/size(err,2);
    case 'avg'
        err=sum(err)/size(err,2);
    case 'med'
        err=median(err);
    otherwise
        error([mth ' not implemented']);
end
end

