function  [logM] = Expectation_v3(deltaX,vacf,PF,trackInfo,stateseq,transmat)

% get track population parameters
numTracks = trackInfo.numberOfTracks;
lambda = trackInfo.lambda;
splitLength = trackInfo.splitLength;
dim = trackInfo.dimensions;
[numStates,numFeatures] = size(vacf);
head=[1;diff(trackInfo.splitIndex)]; %1 means the begining of a track

% calculate analytical covariance terms for each unique track length
invC = cell(1,numStates);
logdetC = zeros(1,numStates);

for k = 1:numStates
    C = toeplitz([vacf(k,:),zeros(1,splitLength-numFeatures-1)]);
    C2 = ShrinkCov(C,lambda);
    [logdetC(k),invC{k}]= LogDeterminant(C2);
end
logM=0;

% calculate normalized posterior probability and log-likelihood
logP = zeros(numTracks,dim);
%Take care of the constant
logP= logP - (splitLength-1)/2*log(2*pi);
for i = 1:numTracks
    state=stateseq(i);
    if head(i)~=0  %every begining of the track, calculate the first term account for transition to itself by
        %choosing the corrsponding poupulation fractio
        P=PF(state);
        if P==0
            P=0.00001;
            disp('pk=0');
        end
        logP(i,:) = logP(i,:)-(.5*diag((deltaX{i}'*invC{state}*deltaX{i})))' + log(P) - .5*logdetC(state);
        
    else
        
        a=stateseq(i-1);
        b=stateseq(i);              %account for state transition
        transP=transmat(a,b);
        
        %To aviod log(0)
        if transP==0
            transP=exp(min(logP(:)));
            disp('transp=0');
        end
        logP(i,:) = logP(i,:)-(.5*diag((deltaX{i}'*invC{state}*deltaX{i})))'+log(transP) - .5*logdetC(state);
    end
end
logM=sum(logP(:));
end
%account for two dimensionss


