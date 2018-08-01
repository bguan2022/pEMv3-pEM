
function  results = pEMv3_SPT(deltaX,trackInfo,params)
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
AccBIC=[];
% BIC Model Selection Loop
state = struct([]);

for numStates = minStates:maxStates
    disp('-------------------------------------------------------');
    disp([num2str(numStates) ' state model']);
    
    % random initialization
    [vacf0,P0] = RandomInitialization(numStates,trackInfo.vacf_exp,2);
    if params.verbose == 1
        disp(['Initial Guess: ' num2str((vacf0(:,1) + 2*vacf0(:,2))'/2/dt)]);
    end
    
    %run rEM
    [baseVacf,baseP,posteriorProb,logLmax,remTrial] = rEM(deltaX,vacf0,P0,params,trackInfo);
    if params.verbose == 1
        disp(['rEM log-likelihood: ' num2str(logLmax)]);
    end
    
    %run pEM
    [baseVacf,baseP,posteriorProb,logLmax,pemTrial] = pEM(deltaX,baseVacf,baseP,params,trackInfo);
    if params.verbose == 1
        disp(['pEM log-likelihood: ' num2str(logLmax)]);
    end
    
    
    disp(['pEM log-likelihood: ' num2str(logLmax)]);
    % run HHM FH
    hmmparams.pik = baseP;
    hmmparams.vacf = baseVacf;
    hmmparams.transProb = .99;
    hmmparams.maxiter = 500;
    hmmparams.tolerance = 1e-3;
    hmmparams.verbose = 1;
    hmmmodel = HMMGMM_v4(deltaX, trackInfo, hmmparams);
    disp(baseP);
    %calculate the state sequence %gamma is the posterior prob, mark the seq by choosing the higher prob
    stateseq = zeros(length(deltaX),1);
    for i = 1:length(deltaX)
        state1 =find(hmmmodel.gammank(i,:) == max(hmmmodel.gammank(i,:)));
        if length(state1)>1
            stateseq(i)=state1(1);
            disp('degnerate');
        else
            stateseq(i)=state1;
        end
        %if found two states of the same prob, pick the first one.
    end
    PF=[];
    est_stateSeq=stateseq;
    l=length(est_stateSeq);
    for k=1:numStates
        pk=sum(est_stateSeq(:) == k)/l;
        PF=[PF pk];
    end
    disp(PF);
    
    HMMlogLmax = hmmmodel.logL;
    bestVacf = hmmmodel.vacf;
    bestP = hmmmodel.PF;
    postProb = hmmmodel.b;
    transmat = hmmmodel.a;

    disp(bestP);
    disp(transmat);
    %calculate logL here after HMM
    [gamma,logL] = Expectation(deltaX,bestVacf,PF,trackInfo);
    newlogL = Expectation_v(deltaX,bestVacf,PF,trackInfo,stateseq,transmat);
    
    % calculate BIC
    nparams = (numStates+1)*(numStates - 1) + numStates*trackInfo.numFeatures; % allows for HMM transition rates
    BIC = newlogL - nparams/2*log(numData);
    AccBIC=[AccBIC BIC];
        disp(['old log L: ' num2str(logL)]);
    disp(['new log L and BIC: ' num2str(newlogL) '  and ' num2str(BIC)]);
    
    %     if BIC(numStates) < max(BIC(1:numStates)) %Since BIC is now negative with HMM,
    %         %max will be nan unless you only look at the BIC you've filled in
    %         break;
    %     end
    
    state(numStates).BIC = BIC; % store results
    state(numStates).numberOfStates = numStates;
    state(numStates).newlogL = newlogL;
    state(numStates).logL = logL;
    state(numStates).vacf = bestVacf;
    state(numStates).bestP = bestP;
    state(numStates).PF = PF;
    state(numStates).posteriorProb = postProb;
    state(numStates).transmat = transmat;
    state(numStates).stateseq = stateseq;
    state(numStates).logLmax = logLmax;
    state(numStates).hmmmodel = hmmmodel;
    
end

% store optimal results
[MAX,numStates] = max(AccBIC);
numStates=numStates+minStates-1;
results.params = params;
results.BIC=AccBIC;
results.state=state;
results.trackInfo = trackInfo;
results.optimalSize = numStates;
results.optimalVacf = state(numStates).vacf;
results.optimalPF = state(numStates).PF;
results.stateseq = state(numStates).stateseq;


toc
end