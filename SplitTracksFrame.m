function [splitX,splitIndex] = SplitTracksFrame(X,splitLength)
%--------------------------------------------------------------------------
% This function splits the tracks into bins of equal lengths, also gives
% frame number corresponding to x y position.  Any remainder
% is not included. The index where the tracks came from is saved in
% splitIndex.
% Feed in X=Xraw, splitlength = 79
%
%
%
% Code written by: 
%       Peter Koo, modified by Mary Lou Bailey on June 3, 2017
%       Yale University, Department of Physics, New Haven, CT, 06511  
%--------------------------------------------------------------------------


%NEED XRAW TO CONTAIN FRAME NUMBERS
n = 1;
counter = 0;
splitX = {};
splitIndex = [];
for i = 1:length(X) %length(X) = 16  
    N = length(X{i}); %length(X{1}) = 90=N
    if N >= splitLength %90>79
        counter = counter + 1; %counter=1
        k = 1;
        status = 1;
        while status == 1
            range = k*splitLength-splitLength+1:k*splitLength;

            if range(end) <= N    
                splitX{n} = X{i}(range,:);
                splitIndex = [splitIndex; counter];
                k = k + 1;
                n = n + 1;
            else
                status = 0;
            end
        end      
    end
end
