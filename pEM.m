function [baseVacf,baseP,posteriorProb,logLmax,trial] = pEM(deltaX,baseVacf,baseP,params,trackInfo)
%--------------------------------------------------------------------------
% This function runs the EM algorithm for a set of perturbation trials, i.e.
% bootstrap set of particle tracks. If the new set of parameters yields a
% higher likelihood with the original data, then the parameters are
% updated.
% 
% Code written by: 
%       Peter Koo
%       Yale University, Department of Physis, New Haven, CT, 06511  
%--------------------------------------------------------------------------


% calculate current likelihood and gamma
[basegamma,logLmax] = Expectation(deltaX,baseVacf,baseP,trackInfo);

% begin perturbation EM to protein trajectories
trial = struct([]);
disp("Iteration in process, please wait...");
for j = 1:params.numPerturbation
%      disp(['EM trial #' num2str(j)]);
    % display current iteration
    
    % perturb likelihood surface with bootstrap resampled data
    trackIndex = ceil(rand(trackInfo.numberOfTracks,1)*trackInfo.numberOfTracks);
    bootStrapDeltaX = cell(trackInfo.numberOfTracks,1);
    vacf_exp = zeros(length(trackIndex),trackInfo.numFeatures,trackInfo.dimensions);
    for i = 1:length(trackIndex)
        bootStrapDeltaX{i} = deltaX{trackIndex(i)};
        vacf_exp(i,:,:) = trackInfo.vacf_exp(trackIndex(i),:,:);
    end

    % calculate new track properties of resampled data
    bsTrackInfo = trackInfo;
    bsTrackInfo.vacf_exp = vacf_exp;
    [vacf_est,P_est,gamma_est,logL] = EM(bootStrapDeltaX,baseVacf,baseP,params,bsTrackInfo);    

    % store results from each trial
    trial(i).vacf_est = vacf_est;
    trial(i).P_est = P_est;
    trial(i).gamma_est = gamma_est;
    trial(i).L = logL;

    % calculate log-likelihood on original data set with bootstrap parameters
    [gamma,Lnew] = Expectation(deltaX,vacf_est,P_est,trackInfo);

    % if bootstrap parameters yield higher likelihood, then employ a full EM with new parameters
    if Lnew > logLmax
        disp('Higher likelihood found');
        [vacf_est,P_est,gamma_est,logL] = EM(deltaX,vacf_est,P_est,params,trackInfo);    
        baseVacf = vacf_est;
        baseP = P_est;
        basegamma = gamma_est;
        logLmax = logL(end);

        trial(i).vacf_est2 = vacf_est;
        trial(i).P_est2 = P_est;
        trial(i).gamma_est2 = gamma_est;
        trial(i).logL2 = logL;
    end    
end

% calculate the posterior probability of each diffusive state
posteriorProb = sum(basegamma,3)/trackInfo.dimensions;




