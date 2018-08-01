function [splitX splitDeltaX splitIndex trackInfo] = SplitTracks2(X,trackInfo)

splitLength = trackInfo.splitLength;
numFeatures = trackInfo.numFeatures;
drift = trackInfo.drift;
rotate = trackInfo.rotate;
trackInfo.trackLength = splitLength - 1;

n = 1;
splitX = {};
splitDeltaX = {};
splitIndex = [];
for i = 1:length(X)  
    
    N = length(X{i});
    k = 1;
    status = 1;
    while status == 1
        range = k*splitLength-splitLength+2-k:k*splitLength+1-k;
        
        if range(end) <= N

            % rotate positions to principal axis
            pos = X{i}(range,:);            
            if rotate == 1
                C = cov(pos);
                [vec val] = eig(C);
                pos2 = (vec*pos')';            
            else
                pos2 = pos;
            end
            
            % calculate properties
            deltaX = diff(pos2);
            [D dim] = size(deltaX);

            splitX{n} = pos2;
            splitDeltaX{n} = deltaX;
            splitIndex = [splitIndex; i];
              
            k = k + 1;
            n = n + 1;
        else
            status = 0;
        end
    end        
end

trackInfo.numberOfTracks = length(splitX);
trackInfo.dimensions = size(splitX{1},2);    % particle track dimensions
trackInfo.splitIndex = splitIndex;

% covariance properties
trackInfo.vacf = CovarianceProperties(splitX,numFeatures,drift);





