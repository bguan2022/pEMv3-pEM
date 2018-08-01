
function  results = pEMv3_SPT_HMMonly(deltaX,trackInfo,params)
%--------------------------------------------------------------------------
% This function runs pEMv2 on a set of particle tracks, X: rEM and then pEM
% to uncover the number of diffusive states and their properties, i.e.
% covariance structure and population fractions. This version assumes all
% tracks share the same length.  
%
% Code written by: 
%       Peter Koo
%       Yale University, Department of Physis, New Haven, CT, 06511  
%--------------------------------------------------------------------------
tic
minStates = params.minStates;
maxStates = params.maxStates;
dt = trackInfo.dt;
numData = length(deltaX)*length(deltaX{1}); 


% BIC Model Selection Loop
state = struct([]);
BIC = zeros(params.maxStates,1); 

for numStates = minStates:maxStates
    disp('-------------------------------------------------------');
    disp([num2str(numStates) ' state model']);

    % random initialization
    [vacf0,P0] = RandomInitialization(numStates,trackInfo.vacf_exp,1);
    if params.verbose == 1
        disp(['Initial Guess: ' num2str((vacf0(:,1) + 2*vacf0(:,2))'/2/dt)]);
    end
    
    %run rEM
    %[baseVacf,baseP,posteriorProb,logLmax,remTrial] = rEM(deltaX,vacf0,P0,params,trackInfo);    
    %if params.verbose == 1
    %    disp(['rEM log-likelihood: ' num2str(logLmax)]);
    %end
    
    %run pEM
    %[baseVacf,baseP,posteriorProb,logLmax,pemTrial] = pEM(deltaX,baseVacf,baseP,params,trackInfo);
    %if params.verbose == 1
    %   disp(['pEM log-likelihood: ' num2str(logLmax)]);
    %end
    
    % run HHM FH
    hmmparams.pik = P0;
    hmmparams.vacfk = Vacf0;
    hmmparams.transProb = .99;
    hmmparams.maxiter = 500;
    hmmparams.tolerance = 1e-7;
    hmmparams.verbose = 1;
    

    hmmmodel = HMMGMM(deltaX, trackInfo, hmmparams);
    
    stateseq = zeros(length(deltaX),1);
    for i = 1:length(deltaX)
        stateseq(i) = find(hmmmodel.gammank(i,:) == max(hmmmodel.gammank(i,:)));
    end
    logLmax = hmmmodel.logL;
    bestVacf = hmmmodel.vacf;
    bestP = hmmmodel.p;
    postProb = hmmmodel.b;
    transmat = hmmmodel.a;
    
    
    % calculate BIC
    % nparams = numStates + numStates*trackInfo.numFeatures; original
    nparams = (numStates+1)*(numStates - 1) + numStates*trackInfo.numFeatures; % allows for HMM transition rates
   %FH numStates-1 since each column of the transition matrix (and the initial population fraction) must add up to 1, 
   %one parameter of a transition matrix column is not an actual parameter (same for the initial population fraction)
   %numStates+1 number of columns in the transition matrix plus the
   %initial population distribution
   
    % BIC(numStates) = logLmax/2 - nparams/2*log(numData); original
    BIC(numStates) = logLmax - nparams/2*log(numData);
    disp(['BIC: ' num2str(BIC(numStates))]);

    % store results
    state(numStates).numberOfStates = numStates;
    state(numStates).BIC = BIC(numStates);
    state(numStates).logL = logLmax;
    state(numStates).vacf = bestVacf';
    state(numStates).P = bestP; 
    state(numStates).posteriorProb = postProb;
    state(numStates).remTrial = remTrial;
    state(numStates).pemTrial = pemTrial;
    state(numStates).transmat = transmat;
    state(numStates).stateseq = stateseq;
    
    if BIC(numStates) < max(BIC(1:numStates)) %Since BIC is now negative with HMM,
        %max will be 0 unless you only look at the BIC you've filled in
    
      break;
    end
end

% store optimal results
[MAX,numStates] = max(BIC(1:numStates));
results.params = params;
results.trackInfo = trackInfo;
results.state = state;
results.optimalSize = numStates;
results.optimalVacf = state(numStates).vacf;
results.optimalP = state(numStates).P;
results.optimalL = state(numStates).logL;
results.posteriorProb = state(numStates).posteriorProb;
results.BIC = BIC;
results.transmat = state(numStates).transmat;
results.stateseq = state(numStates).stateseq;
toc
