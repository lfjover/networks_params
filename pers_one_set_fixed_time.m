function pers_one_set_fixed_time(paramsFile,prevResFile,outDir,iSet)

runName = paramsFile(max(strfind(paramsFile,'/'))+1:end);
runName =runName(min(strfind(runName,'_'))+1:end);
runName = strcat(runName,'_',int2str(iSet));
load(paramsFile)
load(prevResFile)

% matrices = matrices(:,:,50:50:100);
nM = size(matrices,3);

survThres = 10^(-10);
[r,phi,beta,m,x0] = params{iSet,:};

surv = zeros(1,nM);

tic
for iM = 1:nM
    tfinal = 2*stopTimeM(iSet,iM);
    stopTime = tfinal;
    M = matrices(:,:,iM);
    [T,X] = predator_prey_integrator(M,r,a,K,phi,beta,m,x0,tfinal);
    surv(iM) = sum(X(end,:)>=survThres);
    stopTime(iM) = T(end);
%     ['Set: ' int2str(iSet) '. Matrix:' int2str(iM)]
end
runTime = toc;

outFile = strcat(outDir,'result_',runName);
save(outFile,'surv','runTime','stopTime')