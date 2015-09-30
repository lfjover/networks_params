%% Figure 2a: biodiversity vs nestedness ruled-based framework only fm100
close all
clear all

runName = 'rb_fm100_v2';

figure(1);clf;
setfigure(10,6,5,5)
fs = 18;
ms = 6;
run= 1;

par = ['data/params_' runName];
res = ['data/result_' runName];
load(par);load(res)


plot(sortNest(2:end),mean(survM(:,2:end)./20),'ok','markerfacecolor',[0.6 0.6 0.6],...
'markersize',ms)

ylim([0 1])
xlim([sortNest(2) 1])
xlabel('Nestedness (NODF)','interpreter','latex','fontsize',fs)
ylabel('Biodiversity', 'interpreter', 'latex', 'fontsize',fs)


%% Figure 2b: biodiverist vs nestedness ruled-based framework only fm1
runName = 'rb_fm1_v2';

figure(2);clf;setfigure(10,6,5,5)
fs = 18;
ms = 6;
run= 2;

par = ['data/params_' runName];
res = ['data/result_' runName];
load(par);load(res)

% plot persistence
plot(sortNest(2:end),mean(survM(:,2:end)./20),'ok','markerfacecolor',[0.6 0.6 0.6],...
'markersize',ms)
% hold on
% plot(sortNest(focalM),mean(survM(:,focalM)./20),'ok',...
%     'markersize',ms,'markerfacecolor','r')
% hold off
ylim([0 1])
xlim([sortNest(2) 1])
xlabel('Nestedness (NODF)','interpreter','latex','fontsize',fs)
ylabel('Biodiversity', 'interpreter', 'latex', 'fontsize',fs)


%% Figure 3: perturbation-based framework one plot for each trend, single set
runC = {'single_set_fm1_v5','single_set_fm28','single_set_fm100'};


% set fontsize and markersize
fs = 14;
ms = 6;

% set axis properties
figure(3)
clf
px = 18;
py = 6.5;
setfigure(px,py,5,10);
intD = 0.7;
mx = 1.3;
my = 1.2;
axW = (px -2*intD-2*mx)/3;
matH = 1;
matD = 0.2;
axH = py-2*intD-matH-matD;


axes('units','centimeters','position',[ mx my axW axH])
run = 1;
par = ['data/params_' runC{run}]; %parameters (needed for matrices)
res = ['data/result_' runC{run}]; % results
load(par);load(res)
nest = sortNest; %vector with nestedness value
nest(focalM) = []; % erase the nest of focal net to plot as square instead
pers = surv./20;  %biodiviersity (persistence)
pers(focalM) = [];

% fit line
y = surv./20;
p = polyfit(sortNest,y,1);
slope = p(1);
yfit = polyval(p,sortNest);
yresid = y - yfit;
SSresid =  sum(yresid.^2);
SStotal = (length(y)-1)*var(y);
rsq =  (1-SSresid/SStotal);
% calculate p-value
n = length(y);
stdError = sqrt(SSresid/(sum((sortNest - mean(sortNest)).^2)*(n-2)));
tstat = slope/stdError;
pValue_fm1 = 2*tcdf(abs(tstat),n-2,'upper')



plot(nest,pers,'ok','markerfacecolor',[1 .3 .3],...
'markersize',ms)
hold on
plot(sortNest(focalM),surv(focalM)./20,'sk',...
    'markersize',ms,'markerfacecolor',[1 .3 .3])
plot(sortNest,yfit,'-k','linewidth',2)
hold off
ylim([0 1])
xlim([min(sortNest) 1])
xlabel('Nestedness','interpreter','latex','fontsize',fs)
ylabel('Biodiversity', 'interpreter', 'latex', 'fontsize',fs)

% set axis to show fit properties
xtext = 0.79;
ysep = 0.12;
ytop = 0.93;
fsText = 9;
text(xtext,ytop,['\alpha = ' sprintf('%.2f',slope)],...
    'fontsize',fsText)
text(xtext,ytop-ysep,['R^2 = ' sprintf('%.2f',rsq)],...
    'fontsize',fsText)  
% plot matrix
axes('units','centimeters','position',[mx my+axH+matD matH matH])
PlotWeb(matrices(:,:,1))


axes('units','centimeters','position',[ mx+intD+axW my axW axH])
run =2;
par = ['data/params_' runC{run}];
res = ['data/result_' runC{run}];
load(par);load(res)
nest = sortNest;
nest(focalM) = [];
pers = surv./20;
pers(focalM) = [];
y = surv./20;
p = polyfit(sortNest,y,1);
slope = p(1);
yfit = polyval(p,sortNest);
yresid = y - yfit;
SSresid =  sum(yresid.^2);
SStotal = (length(y)-1)*var(y);
rsq =  (1-SSresid/SStotal);
n = length(y);
stdError = sqrt(SSresid/(sum((sortNest - mean(sortNest)).^2)*(n-2)));
tstat = slope/stdError;
pValue_fm28 = 2*tcdf(abs(tstat),n-2,'upper')
plot(nest,pers,'ok','markerfacecolor',[.6 1 .6],...
'markersize',ms)
hold on
plot(sortNest(focalM),surv(focalM)./20,'sk',...
    'markersize',ms,'markerfacecolor',[.6 1 .6])
plot(sortNest,yfit,'-k','linewidth',2)
hold off
set(gca,'ytick',[])
ylim([0 1])
xlim([min(sortNest) 1])
xlabel('Nestedness','interpreter','latex','fontsize',fs)

%fitline
xtext = 0.79;
ysep = 0.12;
ytop = 0.93;
fsText = 9;
text(xtext,ytop,['\alpha = ' sprintf('%.2f',slope)],...
    'fontsize',fsText)
text(xtext,ytop-ysep,['R^2 = ' sprintf('%.2f',rsq)],...
    'fontsize',fsText)  
xPosMat = mx+axW+intD+axW*(sortNest(28)-min(sortNest))/(1-min(sortNest))-matH/2;
axes('units','centimeters','position',...
    [xPosMat my+axH+matD matH matH])
PlotWeb(matrices(:,:,28))

axes('units','centimeters','position',[ mx+2*intD+2*axW my axW axH])
run = 3;
par = ['data/params_' runC{run}];
res = ['data/result_' runC{run}];
load(par);load(res)
nest = sortNest;
nest(focalM) = [];
pers = surv./20;
pers(focalM) = [];
y = surv./20;
p = polyfit(sortNest,y,1);
slope = p(1);
yfit = polyval(p,sortNest);
yresid = y - yfit;
SSresid =  sum(yresid.^2);
SStotal = (length(y)-1)*var(y);
rsq =  (1-SSresid/SStotal);

%p-value
n = length(y);
stdError = sqrt(SSresid/(sum((sortNest - mean(sortNest)).^2)*(n-2)));
tstat = slope/stdError;
pValue_fm100 = 2*tcdf(abs(tstat),n-2,'upper')
plot(nest,pers,'ok','markerfacecolor',[.2 .2 1],...
'markersize',ms)
hold on
plot(sortNest(focalM),surv(focalM)./20,'sk',...
    'markersize',ms,'markerfacecolor',[.2 .2 1])
plot(sortNest,yfit,'-k','linewidth',2)
hold off
set(gca,'ytick',[])
ylim([0 1])
xlim([min(sortNest) 1])
xlabel('Nestedness','interpreter','latex','fontsize',fs)

%fitline
xtext = 0.79;
ysep = 0.12;
ytop = 0.2;
fsText = 9;
text(xtext,ytop,['\alpha = ' sprintf('%.2f',slope)],...
    'fontsize',fsText)
text(xtext,ytop-ysep,['R^2 = ' sprintf('%.2f',rsq)],...
    'fontsize',fsText)  

%matrix
axes('units','centimeters','position',...
    [mx+3*axW+2*intD-matH my+axH+matD matH matH])
PlotWeb(matrices(:,:,100))

  

%% Figure 4a-f: bio vs nest for all deltas for low-nestedness network (fm1)
focalM =1;

par = ['data/params_fm' num2str(focalM) '_delta1_v5'];
load(par)

% create 3d matrix with all the runs;
survM3D = [];
for ii =1 :length(deltaV)
     survFile = ['data/result_fm' num2str(focalM) '_delta' num2str(ii) '_v5'];
    load(survFile) %loads file with survM matrix 
    survM3D = cat(3,survM3D,survM);
end

% fit slope for mean pers for every value of delta
slopes =  [];
rsqV = [];
rSpear = [];
rPear= [];
pValues = [];

figure(4)
setfigure(20,20,60,5)
fs = 18;
ms = 6;
figLetter = {'a','b','c','d','e','f'};
ii = 1;

for dd = 1:size(survM3D,3)
    survM = survM3D(:,:,dd);
    y = mean(survM)./(nH+nV);
    p = polyfit(sortNest,y,1);
    slopes = [slopes p(1)];
    yfit = polyval(p,sortNest);
    yresid = y - yfit;
    SSresid =  sum(yresid.^2);
    SStotal = (length(y)-1)*var(y);
    rsqV = [rsqV (1-SSresid/SStotal)];
    rSpear = [ rSpear corr(sortNest',y','type','Pearson')^2];
    
    %p-value
    n = length(y);
    stdError = sqrt(SSresid/(sum((sortNest - mean(sortNest)).^2)*(n-2)));
    tstat = slopes(end)/stdError;
    pValues = [pValues 2*tcdf(abs(tstat),n-2,'upper')];
       
    subplot(3,2,dd)
    plot(sortNest,y,'ok','markerfacecolor',[1 .3 .3],...
        'markersize',ms)
%     title(['\delta = ' num2str(deltaV(dd)) ' r_S = ' ...
%         num2str(rSpear(end)) ' r_P = ' num2str(rPear(end))...
%         ' m = ' num2str(slopes(end))])
    hold on  
    plot(sortNest,yfit,'-k','linewidth',2)
    hold off
    if dd >4
    xlabel('Nestedness (NODF)','interpreter','latex','fontsize',fs)
    end
    if mod(dd,2) ==1
    ylabel('Biodiversity', 'interpreter', 'latex', 'fontsize',fs)
    end
    ylim([0 1])
    xlim([min(sortNest) max(sortNest)])
    xtext = 0.8;
    ysep = 0.12;
    ytop = 0.93;
    fsText = 12;
    text(xtext,ytop,['\delta = ' num2str(deltaV(dd))],'fontsize',fsText)
    text(xtext,ytop-ysep,['\alpha = ' sprintf('%.2f',slopes(end))],...
        'fontsize',fsText)
    text(xtext,ytop-2*ysep,['R^2 = ' sprintf('%.2f',rsqV(end))],...
        'fontsize',fsText) 
    text(min(sortNest),1.12,figLetter{ii},'fontsize',20)
    ii = ii+1;
end


%% Figure 4g-l: bio vs nest for all deltas for  perfectly nested net (fm100)
focalM =100; %focal matrix
par = ['data/params_fm' num2str(focalM) '_delta1_v4'];
load(par)


% create 3d matrix with all the runs ( all values of delta);
survM3D = [];
for ii =1 :length(deltaV)
     survFile = ['data/result_fm' num2str(focalM) '_delta' num2str(ii) '_v4'];
    load(survFile) %loads file with survM matrix 
    survM3D = cat(3,survM3D,survM);
end


slopes =  [];
rsqV = [];
rSpear = [];
rPear= [];
pValues = [];

figure(5)
setfigure(20,20,10,5)
fs = 18;
ms = 6;

figLetter = {'g','h','i','j','k','l'};
ii = 1;
for dd = 1:size(survM3D,3)
    % fit slope for mean biodiversity for every value of delta
    survM = survM3D(:,:,dd);
    y = mean(survM)./(nH+nV);
    p = polyfit(sortNest,y,1);
    slopes = [slopes p(1)];
    yfit = polyval(p,sortNest);
    yresid = y - yfit;
    SSresid =  sum(yresid.^2);
    SStotal = (length(y)-1)*var(y);
    rsqV = [rsqV (1-SSresid/SStotal)];
    rSpear = [ rSpear corr(sortNest',y','type','Pearson')^2];
    
    %p-value
    n = length(y);
    stdError = sqrt(SSresid/(sum((sortNest - mean(sortNest)).^2)*(n-2)));
    tstat = slopes(end)/stdError;
    pValues = [pValues 2*tcdf(abs(tstat),n-2,'upper')];
       
    subplot(3,2,dd)
    plot(sortNest,y,'ok','markerfacecolor',[.2 .2 1],...
        'markersize',ms)
%     title(['\delta = ' num2str(deltaV(dd)) ' r_S = ' ...
%         num2str(rSpear(end)) ' r_P = ' num2str(rPear(end))...
%         ' m = ' num2str(slopes(end))])
    hold on  
    plot(sortNest,yfit,'-k','linewidth',2)
    hold off
    if dd >4
    xlabel('Nestedness (NODF)','interpreter','latex','fontsize',fs)
    end
%     if mod(dd,2) ==1
%     ylabel('Biodiversity', 'interpreter', 'latex', 'fontsize',fs)
%     end
    ylim([0 1])
    xlim([min(sortNest) max(sortNest)])
    xtext = 0.82;
    ysep = 0.12;
    ytop = 0.3;
    fsText = 12;
    text(xtext,ytop,['\delta = ' num2str(deltaV(dd))],'fontsize',fsText)
    text(xtext,ytop-ysep,['\alpha = ' sprintf('%.2f',slopes(end))],...
        'fontsize',fsText)
    text(xtext,ytop-2*ysep,['R^2 = ' sprintf('%.2f',rsqV(end))],...
        'fontsize',fsText)  
    text(min(sortNest),1.12,figLetter{ii},'fontsize',20)
    ii = ii+1;
end

%% Figure 6: Focal networks in perturbation-based framework
load('matrices/matrices_10by10_100_invertible.mat')
fs =18;
figure(6)
setfigure(20,7,10,10)
subplot(1,3,3)
M = matrices(:,:,end);
PlotWeb(M)
title('virus','fontsize',fs,'interpreter','latex')
ylabel('bacteria','fontsize',fs,'interpreter','latex')
xlabel(['NODF = ' num2str(sortNest(end),'%0.2f')],...
        'fontsize',fs-2,'interpreter','latex')
subplot(1,3,2)
num =28;
M = matrices(:,:,num);
PlotWeb(M)
title('virus','fontsize',fs,'interpreter','latex')
ylabel('bacteria','fontsize',fs,'interpreter','latex')
xlabel(['NODF = ' num2str(sortNest(num),'%0.2f')],...
        'fontsize',fs-2,'interpreter','latex')
subplot(1,3,1)
M = matrices(:,:,1);
PlotWeb(M)
title('virus','fontsize',fs,'interpreter','latex')
ylabel('bacteria','fontsize',fs,'interpreter','latex')
xlabel(['NODF = ' num2str(sortNest(1),'%0.2f')],...
        'fontsize',fs-2,'interpreter','latex')

%% Figure 7: focal networks in ruled-based framework 
%  low nestedness and perfectly nested networks
fs =20;
figure(7)
setfigure(13,5,10,10)
load('matrices/matrices_10by10_100_special.mat')
subplot(1,2,1)
M = matrices(:,:,1);
PlotWeb(M)
title('virus','fontsize',fs,'interpreter','latex')
ylabel('bacteria','fontsize',fs,'interpreter','latex')
text(-1.5,12.5,'a','fontsize',fs)

subplot(1,2,2)
M = matrices(:,:,end);
PlotWeb(M)
title('virus','fontsize',fs,'interpreter','latex')
ylabel('bacteria','fontsize',fs,'interpreter','latex')
text(-1.5,12.5,'b','fontsize',fs)

%% Fgure 8; doubling time
figure(8);clf;
setfigure(20,15,5,5)
fs = 14;
ms = 6;

par = 'data/params_rb_fm100_v2';
res = 'data/result_rb_fm100_ft';
load(par);load(res)
sm = survM;

% plot persistence
subplot(2,2,2)
plot(sortNest,mean(survM./20),'ok','markerfacecolor',[0.6 0.6 0.6],...
'markersize',ms)
% hold on
% plot(sortNest(focalM),mean(survM(:,focalM)./20),'ok',...
%     'markersize',ms,'markerfacecolor','r')
% hold off
ylim([0 1])
xlabel('Nestedness (NODF)','interpreter','latex','fontsize',fs)
ylabel('Biodiversity', 'interpreter', 'latex', 'fontsize',fs)
%title('Double time')


res = 'data/result_rb_fm100_v2';
load(par);load(res)

% plot persistence
subplot(2,2,1)
plot(sortNest,mean(survM./20),'ok','markerfacecolor',[0.6 0.6 0.6],...
'markersize',ms)
% hold on
% plot(sortNest(focalM),mean(survM(:,focalM)./20),'ok',...
%     'markersize',ms,'markerfacecolor','r')
% hold off
ylim([0 1])
xlabel('Nestedness (NODF)','interpreter','latex','fontsize',fs)
ylabel('Biodiversity', 'interpreter', 'latex', 'fontsize',fs)
%title('Stopping time heuristics')

subplot(2,1,2)
hist((mean(sm./20)-mean(survM./20))./mean(survM./20)*100)
xlabel('Percentage change in average biodiversity')
