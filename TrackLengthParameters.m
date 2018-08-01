function [trackLength,uniqueLength] = TrackLengthParameters(deltaX)
%-------------------------------------------------------------------------- 
% This function finds the number of displacements for each
% particle track and also finds unique values of the track length.  This
% helps reduce unnecessary recalculations of the inverse and determinants of 
% particle tracks with the same lengths.
%
% Code written by: 
%       Peter Koo
%       Yale University, Department of Physis, New Haven, CT, 06511  
%-------------------------------------------------------------------------- 

numTracks = length(deltaX);
trackLength = zeros(numTracks,1);
for i = 1:numTracks
	trackLength(i) = size(deltaX{i},1);
end
uniqueLength = unique(trackLength);

