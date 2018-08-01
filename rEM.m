function [vacf,P,posteriorProb,logLmax,trial] = rEM(deltaX,vacf0,P0,params,trackInfo)
%--------------------------------------------------------------------------
% This function runs the EM algorithm with different parameter initializations
% and returns the parameters from the trial which yield the highest 
% likelihood value
%
% Code written by: 
%       Peter Koo
%       Yale University, Department of Physis, New Haven, CT, 06511  
%--------------------------------------------------------------------------

numTrial = params.numReinitialize;

% initialization
[numStates,numFeatures] = size(vacf0);
baseVacf = zeros(numStates,numFeatures,numTrial);
baseP = zeros(numTrial,numStates);
basegamma = cell(numTrial);
logLmax = zeros(numTrial,1);
trial = struct([]);

% reinitialization EM trials
for i = 1:numTrial
    disp(['EM trial #' num2str(i)]);

    % run em algorithm
    [vacf_est,P_est,gamma_est,logL] = EM(deltaX,vacf0,P0,params,trackInfo);
    try
        baseVacf(:,:,i) = vacf_est;
        baseP(i,:) = P_est;
        basegamma{i} = gamma_est;
        logLmax(i) = logL(end);
    end
    
    % store results from each trial
    trial(i).vacf_est = vacf_est;
    trial(i).P_est = P_est;
    trial(i).gamma_est = gamma_est;
    trial(i).L = logL;
    
    % reinitialize EM parameters
    [vacf0,P0] = RandomInitialization(numStates,trackInfo.vacf_exp,2);
    
end

% find maximum likelihood across trials
[MAX,index] = max(logLmax);
vacf  = baseVacf(:,:,index);
P  = baseP(index,:);
logLmax = MAX;

% calculate the posterior probability of each diffusive state
gamma = basegamma{index};
posteriorProb = sum(gamma,3)/trackInfo.dimensions;




