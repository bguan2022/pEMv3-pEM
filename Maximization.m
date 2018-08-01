function [vacf,P] = Maximization(gamma,vacf_exp)
%--------------------------------------------------------------------------
% This function calculates the maximization step: the covariance elements and
% the population fractions given the posterior probabilities.
% 
% Code written by: 
%       Peter Koo
%       Yale University, Department of Physis, New Haven, CT, 06511  
%--------------------------------------------------------------------------

% get track population parameters
[numTracks,numStates,dim] = size(gamma);
numFeatures = size(vacf_exp,2);

% population fraction maximization
Nk = zeros(dim,numStates);
for i = 1:dim
    Nk(i,:) = sum(gamma(:,:,i));
end
meanNk = mean(Nk,1);
P = meanNk/numTracks;


% diffusion coefficient maximization
vacf = zeros(numStates,numFeatures);
for j = 1:dim
    for k = 1:numFeatures
        vacf(:,k) = vacf(:,k) + (sum(gamma(:,:,j).*(vacf_exp(:,k,j)*ones(1,numStates)))./Nk(j,:))';
    end
end
vacf = vacf/dim;

