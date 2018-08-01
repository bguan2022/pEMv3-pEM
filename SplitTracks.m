function [splitX,splitIndex] = SplitTracks(X,splitLength)
%--------------------------------------------------------------------------
% This function splits the tracks into bins of equal lengths.  Any remainder
% is not included. The index where the tracks came from is saved in
% splitIndex.
%
% Code written by: 
%       Peter Koo
%       Yale University, Department of Physis, New Haven, CT, 06511  
%--------------------------------------------------------------------------

n = 1;
counter = 0;
splitX = {};
splitIndex = [];
for i = 1:length(X)  
    N = length(X{i});
    if N >= splitLength
        counter = counter + 1;
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
