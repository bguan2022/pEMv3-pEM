

clear all;
close all;
addpath('/Users/Jguan/Desktop/Yale/pem_v3');
addpath('/Users/Jguan/Desktop/Yale/pEM_v3/pEMv3');
addpath('/Users/Jguan/Desktop/Yale/pEM_v3/HMM');

%% 1. load file
load('/Users/Jguan/Desktop/Yale/pEMv3/simulations_new/data/200tracks/case1.mat');
savename = 'case1.mat';

%% 2. user set parameters
% movie parameters
dt = .058;              % time between steps
dE = .01;              % exposure time

% pEM parameters
minStates =4;          % minimum number of states to explore
maxStates = 5;          % maximum number of states to explore
numReinitialize = 12;    % was 6 -- number of reinitialization trials
numPerturb = 300;        % was 100 -- number of perturbation trials
maxiter = 1000;          % maximum number of iterations within EM trial
convergence = 1e-5;     % 1e-5 convergence criteria for change in log-likelihood
lambda = 0.0001;        % shrinkage factor (useful when numerical issues calculating
% inverse of covariance matrix, labmda = 0.0 for no correction
% lambda = 0.001 for correction)

splitLength=20;        %must be at least numFeatures+1
numFeatures =splitLength-1 ;
minLength = splitLength;

plotMSD=1;
%% load files 
%
file1='/Users/Jguan/Desktop/Yale/pEMv3/Tet/TetR_SMT/Cell_4b_evalSPT-TRACES.txt';
file2='/Users/Jguan/Desktop/Yale/pEMv3/Tet/TetR_SMT/Cell_9a_evalSPT-TRACES.txt';
file3='/Users/Jguan/Desktop/Yale/pEMv3/Tet/TetR_SMT/Cell_9b_evalSPT-TRACES.txt';
file4='/Users/Jguan/Desktop/Yale/pEMv3/Tet/TetR_SMT/Cell_10_evalSPT-TRACES.txt';
file5='/Users/Jguan/Desktop/Yale/pEMv3/Tet/TetR_SMT/Cell_11_evalSPT-TRACES.txt';
file6='/Users/Jguan/Desktop/Yale/pEMv3/Tet/TetR_SMT/Cell_12_evalSPT-TRACES.txt';
file7='/Users/Jguan/Desktop/Yale/pEMv3/Tet/TetR_SMT/Cell_13a_evalSPT-TRACES.txt';
file8='/Users/Jguan/Desktop/Yale/pEMv3/Tet/TetR_SMT/Cell_13b_evalSPT-TRACES.txt';

% % file9='/Users/Jguan/Desktop/Yale/pEMv3/Tet/TetO_Coinjection/Cell_4-TetO_TRACES.txt';
% file10='/Users/Jguan/Desktop/Yale/pEMv3/Tet/TetO_Coinjection/Cell_4-TetO_TRACES.txt';
% file11='/Users/Jguan/Desktop/Yale/pEMv3/Tet/TetO_Coinjection/Cell_4-TetO_TRACES.txt';
% file12='/Users/Jguan/Desktop/Yale/pEMv3/Tet/TetO_Coinjection/Cell_4-TetO_TRACES.txt';
% %
% % file1=('/Users/Jguan/Desktop/Yale/pEMv3/simulations_new/data/715/DeltaX_state1.txt');
% % file2=('/Users/Jguan/Desktop/Yale/pEMv3/simulations_new/data/715/DeltaX_state2.txt');
% % file3=('/Users/Jguan/Desktop/Yale/pEMv3/simulations_new/data/715/DeltaX_state3.txt');
% % file4=('/Users/Jguan/Desktop/Yale/pEMv3/simulations_new/data/715/DeltaX_state4.txt');
% % file5=('/Users/Jguan/Desktop/Yale/pEMv3/simulations_new/data/715/DeltaX_state5.txt');
% % % phena
file = {file1,file2, file3, file4, file5, file6, file7, file8};
% file = {file9,file10, file11, file12};

pixels2um = 1.0;
formatspec='%f  %f      %d      %d      %f      %f      %f      %f';
oldtracks=0;

for j=1:length(file) %loop through each file, change from 1:14 to 1:10?
    fid=fopen(file{j},'r');
    data=textscan(fid,formatspec); %cell containing each column from the text file
    fclose(fid);
    A=data{4}; %get 4th cell - track numbers
    B=unique(A); %get one of each in A
    tracks=length(B); %total number of tracks found in the 2D code
    for i=1:tracks
        index=find(A==B(i));
        X{i+oldtracks}=[pixels2um*data{1}(index),pixels2um*data{2}(index)]; %save x, y positions
        sigMA{i+oldtracks}=data{4}(index); %save sigmas
        trackIndex(i+oldtracks)=j;
    end
    oldtracks= oldtracks+tracks;
end
%cell containing the x and y positions for each track
%% 3. Can only use tracks longer than splitLength
Xraw=X;
X2={};
k = 1;
% Xraw=AllX
for i = 1:length(Xraw) %loop through each track, make sure it's at least splitlength long
    if size(Xraw{i},1) >10
        X2{k} = Xraw{i};
        %         sigMA2{k}=sigMA{i};
        %         trackIndex2(k) = trackIndex(i);
        k = k + 1;
    end
end
% trackIndex=trackIndex2; %list of corresponding track numbers
Xraw = X2;
disp(['Using ' num2str(k-1) ' tracks >= ' num2str(minLength) ' steps long.']);



%% 4. run pEM version 3
% split tracks into equal bin sizes
% clear X;
[X,splitIndex] = SplitTracks(Xraw,splitLength); %feed in x y positions for each track, splitlength,
%returns x y positions

% structure for track info
trackInfo.numberOfTracks = length(X);   % number of tracks
trackInfo.dimensions = size(X{1},2);    % particle track dimensions
trackInfo.numFeatures = numFeatures;    % number of features to retain in covariance matrix
trackInfo.splitLength = splitLength;    % length of each bin
trackInfo.splitIndex = splitIndex;      % index of each track
trackInfo.dt = dt;                      % frame duration
trackInfo.R = 1/6*dE/dt;                % motion blur coefficient
trackInfo.lambda = lambda;              % shrinkage factor
[lengthoftracks,dim]=size(Xraw{1});

% trackInfo.markovStateSeq=markovStateSeq;
% trackInfo.transmat=simParams.A;

% structure for pEM
params.minStates = minStates;               % minimum number of states to try
params.maxStates = maxStates;               % maximum number of states to try
params.numFeatures = numFeatures;           % number of features in covariance elements
params.numReinitialize = numReinitialize;   % number of reinitialization trials
params.numPerturbation = numPerturb;        % number of perturbations trials
params.converged = convergence;             % convergence condition for EM
params.maxiter = maxiter;                   % maximum number of iterations for EM
params.verbose = 0;                         % display progress on command window (0,1)


% calculate the displacements for each particle track
deltaX = cell(trackInfo.numberOfTracks,1);
for i = 1:trackInfo.numberOfTracks
    deltaX{i} = diff(X{i});
end

% calculate relevant properties to enhance compuatational time
[trackInfo.vacf_exp,trackInfo.xbar_exp] = CovarianceProperties(deltaX,numFeatures);

%%
% run pEMv3
results = pEMv3_SPT(deltaX,trackInfo,params);


%%
% display results
optimalSize = results.optimalSize;
optimalVacf = results.state(optimalSize).vacf;
%results.optimalP gives optimal inital pi for each track, so
%take the average of the optimal inital pi's


% optimalVacf = results.state(numStates).vacf;
% optimalP = results.state(numStates).P;
% optimalPF = results.state(numStates).PF;
% optimalL = results.state(numStates).logL;

% est_stateSeq = results.state(numStates).stateseq;
% transmat=results.state(numStates).transmat;

% disp('-------------------------------------------------------');
% disp(['OptimalSize: ' num2str(optimalSize) ' states']);
% for i = 1:optimalSize
%     disp(['Sigma_k(i,i+' num2str(i-1) '): ' num2str(optimalVacf(:,i)') ' um^2']);
% end

%%
% est_stateSeq=est_stateIndex_orig; %pEMv2
est_stateSeq=results.state(optimalSize).stateseq;
% real_stateSeq=cat(1,markovStateSeq{:});
% real_stateSeq=real_stateSeq(1:splitLength:end);

% Rate(real_stateSeq,est_stateSeq);

%% 7. Setup variables
results.X=X;
splitX = results.X;
trackInfo = results.trackInfo;
splitIndex = results.trackInfo.splitIndex;
splitLength = trackInfo.splitLength;
numFeatures = trackInfo.numFeatures;
vacf_exp = trackInfo.vacf_exp;
BIC = results.BIC;
state = results.state; %logL is stored in state
numTracks=trackInfo.numberOfTracks;
numStates = results.optimalSize;
A=results.state(optimalSize).transmat;
%% 8. Print summary of results
disp('-------------------------------------------------------');
disp(['Summary of Results']);
disp(['Est number of states: ' num2str(results.optimalSize) ]);

%% 10. print estimated transition matrix: NEED TO CACULATE TRANSITION MATRIX IN VERSION 3
disp('estimated HMM transition matrix:');
disp(num2str(A));

%% 11. Estimated MSD and  Estimated covariance structure for each state, plots MSD and Vacf

plotMSD=1;
[MSD,Vacf]=MSD_Vacf(X,est_stateSeq,plotMSD,6);

%%
% store results in analytics structure
analytics = struct;
% analytics.MSD = MSD;
% analytics.Vacf = Vacf;
% analytics.BIC = BIC;
% analytics.est_stateSeq = est_stateSeq;

%% 16 Save variables to mat file
save(strcat('/Users/Jguan/Desktop/Yale/pEMv3/simulations_new/data/805/X_',savename),'X','Xraw','A','analytics','results','Xraw','splitIndex','est_stateSeq','BIC')
disp('results saved');

%% 17. Plot tracks
% plotTracks4(Xraw,markovStateSeq,stateSeq); %after pEM
% plotTracks1(Xraw,markovStateSeq); %before pEM
% plotTracks0(X');

%% %PLOT LogL
clear BIC;clear numSatates;
min=1;max=8;



newlogL=[results.state.newlogL];
logL=[results.state.logL];

BIC=[results.BIC];
oldBIC=[];
for numStates = 1:8     % calculate BIC
    logM=logL(numStates);
    nparams = (numStates+1)*(numStates - 1) + numStates*numFeatures; % allows for HMM transition rates

    BIC1 = logM- nparams/2*log(length(results.stateseq)*5);
    oldBIC=[oldBIC BIC1];

end
% % plot BIC

figure; hold on;
plot([min:max],newlogL,'r');
plot([min:max],BIC,'k');

plot([min:max],oldBIC,'b');
plot([min:max],logL,'g');
% plot([min:max],logLmax,'g');

xlabel('Model size (states)','fontsize',20);
ylabel('LogL','fontsize',20);

hold on;
h = zeros(4, 1);
h(1) = plot(NaN,NaN,'-r');
h(2) = plot(NaN,NaN,'-k');
h(3) = plot(NaN,NaN,'-g');
h(4) = plot(NaN,NaN,'-b');

legend(h,{'logL','BIC','oldlogL','oldBIC'},'fontsize',16);
% Lifetime(stateSeq,splitLength);

