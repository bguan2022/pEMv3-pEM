function  [gamma,logL] = Expectation(deltaX,vacf,PF,trackInfo)
%--------------------------------------------------------------------------
% This function calculates the expectation step: the posterior probabilities
% and the log-likelihood given the model parameters
% 
% Code written by: 
%       Peter Koo
%       Yale University, Department of Physis, New Haven, CT, 06511  
%--------------------------------------------------------------------------

% get track population parameters
numTracks = trackInfo.numberOfTracks;
lambda = trackInfo.lambda;
splitLength = trackInfo.splitLength;
dim = trackInfo.dimensions;
[numStates,numFeatures] = size(vacf);


% calculate analytical covariance terms for each unique track length
invC = cell(1,numStates);
logdetC = zeros(1,numStates);
for k = 1:numStates
    C = toeplitz([vacf(k,:),zeros(1,splitLength-numFeatures-1)]);
    C2 = ShrinkCov(C,lambda);
    [logdetC(k),invC{k}]= LogDeterminant(C2);        
end

% calculate normalized posterior probability and log-likelihood
logP = zeros(numTracks,numStates,dim); 
for k = 1:numStates
    for i = 1:numTracks
            logP(i,k,:) = -.5*diag((deltaX{i}'*invC{k}*deltaX{i}));
    end
end

logL = 0; 

gamma = zeros(numTracks,numStates,dim); 
for j = 1:dim 
    logpiN = logP(:,:,j) + ones(numTracks,1)*log(PF) - (splitLength-1)/2*log(2*pi) - .5*ones(numTracks,1)*logdetC;

    % fix numerical issues associated with log-sum-exp
    MAX = max(logpiN,[],2); %find the maximun probability along the second dimension
    logpiN_norm = bsxfun(@minus,logpiN,MAX);
    logSumPiN = MAX + log(sum(exp(logpiN_norm),2));
    index = find(~isfinite(MAX));
    if ~isempty(index)
        logSumPiN(index) = MAX(index);
    end
    logr = logpiN - repmat(logSumPiN,1,numStates);
    gammank = exp(logr);
    gammank(find(gammank(:) < 1e-100)) = 1e-100;
    gamma(:,:,j) = gammank; 
    logL = logL + sum(logSumPiN);
 
end


