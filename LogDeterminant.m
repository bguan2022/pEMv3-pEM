function [logdetC,invC]= LogDeterminant(C)
%-------------------------------------------------------------------------- 
% Summary: This function calculates the log determinant and inverse of the
% covariance matrix with eigenvalues to circumvent numerical issues
% 
% Input:
%       C = covariance matrix
%
% Output:
%       logdetC = log determinant of covariance matrix
%       invC = inverse of covariance matrix
%
% Code written by: 
%       Peter Koo
%       Yale University, Department of Physis, New Haven, CT, 06511  
%-------------------------------------------------------------------------- 

threshold = 1e-10;
 C(isnan(C)) = 0.000001;
 C(isinf(C)) = 0.000001;
[vec val] = eig(C);
val = diag(val);
index = find(val <=threshold);
val(index) = threshold;
logdetC = sum(log(val));
invval = diag(1./val);
invC = vec*invval*vec';

end