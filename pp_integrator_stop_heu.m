function [T,X] = pp_integrator_stop_heu(M,r,a,K,phi,beta,m,x0,timeInc,...
                                   totalTime,survThres,rErrTh)
% PP_INTEGRATOR_STOP_HEU integration of predator_pray dynamics,stops
% if surviving species(densities greater than threshold) are rErrTh close
% to theory or if dynamics go for totalTime. Conditions are checked every 
% timeInc.


rErr = 1;
tfinal = 0;
T = [0];
X = [x0'];
xEqui = 0;

while any(abs(rErr)>rErrTh) && tfinal < totalTime
        
    [t,x] = predator_prey_integrator(M,r,a,K,phi,beta,m,x0,timeInc);
    T = [T;t(2:end)+tfinal];
    X = [X;x(2:end,:)];
    tfinal = tfinal+timeInc;
    x0 = x(end,:)';
    
    iSurv = x(end,:)>survThres;
    Msurv = M(iSurv(1:10),iSurv(11:20));
    [nHsurv,nVsurv] = size(Msurv);
        
    if rank(Msurv)==max([nHsurv nVsurv])
        [rSurv,phiSurv,betaSurv,mSurv,aSurv] = get_surv_params(iSurv,M,r,phi,beta,m);
        xEqui = equilibrium(Msurv,rSurv,aSurv,K,phiSurv,betaSurv,mSurv);
        xSurv = X(floor(length(T)/2):end,iSurv);
        xNum = mean(xSurv,1);
        rErr=(xEqui-xNum')./xEqui;
    end
end

% iSurv = x(end,:)>survThres;
% iV = logical([zeros(1,10) iSurv(11:20)]);
% iH = logical([iSurv(1:10) zeros(1,10)]);
% figure(1);
% subplot(2,1,1);
% semilogy(T(floor(length(T)*8/10):end),X(floor(length(T)*8/10):end,iSurv(1:10)));
% title('host')
% % legend
% % hold on; 
% % if xEqui ~= 0
% % semilogy([T(floor(length(T)*8/10)) T(end)],repmat(xEqui(1:10),1,2),'k-'); 
% % end
% % hold off
% subplot(2,1,2);
% semilogy(T(floor(length(T)*8/10):end),X(floor(length(T)*8/10):end,iV));
% % hold on; 
% % if xEqui ~= 0
% % semilogy([T(floor(length(T)*8/10)) T(end)],repmat(xEqui(11:20),1,2),'k-');
% % end
% % hold off
% figure(2)
% subplot(1,2,1);PlotWebColor(M,find(iSurv(1:10)),find(iSurv(11:20)),'r')
% subplot(1,2,2);PlotWeb(M(iSurv(1:10),iSurv(11:20)))
% figure(3)
% semilogy(T(1:50:end),X(1:50:end,iSurv))
% figure(4)
% semilogy(T(1:100:end),X(1:100:end,:))

function [rSurv,phiSurv,betaSurv,mSurv,aSurv] = get_surv_params(iSurv,M,r,phi,beta,m)
    
    Hsurv = iSurv(1:10);
    Vsurv = iSurv(11:20);
    Msurv = M(Hsurv,Vsurv);
    [nHsurv,~] = size(Msurv);
    
    rSurv = r(Hsurv);
    phiSurv = phi(1:nHsurv,Vsurv);
    betaSurv = beta(1:nHsurv,Vsurv);
    mSurv = m(Vsurv);
    aSurv = ones(nHsurv);