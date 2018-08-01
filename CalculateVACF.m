function vacf = CalculateVACF(C)
%--------------------------------------------------------------------------
% This function calculates the covariance structure of matrix C assuming a
% symmetric toeplitz matrix.
%
% Code written by: 
%       Peter Koo
%       Yale University, Department of Physis, New Haven, CT, 06511  
%--------------------------------------------------------------------------

D = length(C);
subC = zeros(D);
for i = 1:D
    subC(i,1:(D-i+1)) = C(i,i:end);
end

vacf = sum(subC)./(D:-1:1);
