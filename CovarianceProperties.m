function [vacf_exp,xbar_exp] = CovarianceProperties(deltaX,numFeatures)
%--------------------------------------------------------------------------
% This function calculates the covariance structure and mean of each particle 
% track.  the covariance structure is only calculated to the number of
% features and places a zero beyond.
%
% Code written by: 
%       Peter Koo
%       Yale University, Department of Physis, New Haven, CT, 06511  
%--------------------------------------------------------------------------

numTracks = length(deltaX);
dim = size(deltaX{1},2);

% calculate velocity auto-covariance for all particle tracks
vacf_exp = zeros(numTracks,numFeatures,dim); 
for j = 1:dim
    for i = 1:numTracks
        diffX = deltaX{i}(:,j);            
        C = diffX*diffX';
        tmpvacf = CalculateVACF(C);
        vacf_exp(i,:,j) = tmpvacf(1:numFeatures);
    end
end 

% mean
xbar_exp = zeros(numTracks,dim);
for i = 1:numTracks
    xbar_exp(i,:) = mean(deltaX{i});
end


