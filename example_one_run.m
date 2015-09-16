
runName = 'single_set_fm1';
paramsFile = ['data/params_' runName];
outDir = 'data/example/';
iSet = 1;

% calculate biodiversity for 100 matrices
pers_one_set(paramsFile,outDir,iSet)


%% plot biodiversity vs. nestedness
res = [outDir 'result_' runName '_' num2str(iSet)];
load(res) %load biodiversity result
load(paramsFile) %load nestedness values for matrix ensemble

fs = 18;
ms = 6;
plot(sortNest,surv./20,'ok',...
    'markerfacecolor',[0.6 0.6 0.6],'markersize',ms)
ylim([0 1])
xlim([sortNest(1) 1])
xlabel('Nestedness (NODF)','interpreter','latex','fontsize',fs)
ylabel('Biodiversity', 'interpreter', 'latex', 'fontsize',fs)
setfigure(10,8,6,6);print('-dpdf','bio_one_set.pdf')
setfigure(10,8,6,6);print('-dpng','bio_one_set.png')