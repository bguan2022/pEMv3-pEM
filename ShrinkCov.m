function C2 = ShrinkCov(C, lambda)
%--------------------------------------------------------------------------
% This function calculates the Ledoit-Wolf optimal shrinkage estimator for 
% the covariance matrix. This is only necessary when inversion of the
% covariance matrix is problematic. 
% 
% Code written by: 
%       Peter Koo
%       Yale University, Department of Physis, New Haven, CT, 06511  
%--------------------------------------------------------------------------

% shrink C
C2 = (1-lambda)*C + lambda*diag(diag(C));
        
end
