function pers_one_set(paramsFile,outDir,iSet)
% integrates the system for all the matrices (100) in paramsFile and the
% parameter set number iSet in paramsFile
% output: 
% surv: number of species with densities > threshold for each matrix.
% stopTime: stopping time for each simulation
% runTime: total run time of 100 matrices

%create the name of the run (used in output file)
runName = paramsFile(max(strfind(paramsFile,'/'))+1:end);
runName =runName(min(strfind(runName,'_'))+1:end);
runName = strcat(runName,'_',int2str(iSet));

%load the file with model parameters
load(paramsFile)

totalTime = 40000; % total time of simulation if other condition isn't 
                   % reached
timeInc = 500;
survThres = 10^(-10);
rErrTh = 0.1;
nM = size(matrices,3);

surv = zeros(1,nM);
stopTime= zeros(1,nM);

[r,phi,beta,m,x0] = params{iSet,:}; % assign params to variables


% integrate the system
tic
for iM = 1:nM
    M = matrices(:,:,iM); %select one matrix out of 100.
    [T,X] = pp_integrator_stop_heu(M,r,a,K,phi,beta,m,x0,timeInc,...
                           totalTime,survThres,rErrTh);
    surv(iM) = sum(X(end,:)>=survThres); %calculate survivors
    stopTime(iM) = T(end);
    ['Matrix: ' int2str(iM) '/' int2str(nM)]
end
runTime = toc;

outFile = strcat(outDir,'result_',runName);
save(outFile,'surv','stopTime','runTime')