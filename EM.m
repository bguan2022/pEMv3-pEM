function  [vacf_est,P_est,gamma_est,logL] = EM(deltaX,vacf,P,params,trackInfo)
%--------------------------------------------------------------------------
% This function runs the EM algorithm for a set of particle track
% displacements, deltaX, and outputs the converged diffusion model
% parameters for each diffusive state, i.e. vacf and population fraction,
% posterior probability, and log-likelihood value.
% 
% Code written by: 
%       Peter Koo
%       Yale University, Department of Physis, New Haven, CT, 06511  
%--------------------------------------------------------------------------

% parameters
converged = params.converged;
maxiter = params.maxiter;
vacf_exp = trackInfo.vacf_exp;

% expectation step with a initial parameters
[gamma,logL] = Expectation(deltaX,vacf,P,trackInfo);

% expectation-maximization at each annealing temperature
vacf_est = vacf; P_est = P; gamma_est = gamma;
for i = 2:maxiter    
    
    % maximization step
    [vacf,P] = Maximization(gamma,vacf_exp);

    % expectation step    
    [gamma,logL(i,1)] = Expectation(deltaX,vacf,P,trackInfo);
        
    % check for convergence of the log-likelihood
    if logL(i) - logL(i-1) < converged 
        logL(end) = [];
        break;
    else
        % if not converged, then store updated parameters
        vacf_est = vacf;
        P_est = P;
        gamma_est = gamma;
    end
end 






