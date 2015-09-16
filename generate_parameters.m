%% parameters for the rule-based framework and perfectly nested net (fm100)
clear all
load('matrices/matrices_10by10_100_special_inv');
focalM = 100;
name = ['params_rb_fm' int2str(focalM) '_v2'];

nSet = 100;

[nH,nV] = size(matrices(:,:,1));
a = ones(nH);

rMin = 3.61; rMax=43.35;
phiMin = 2.4*10^-7; phiMax = 2.4*10^-6;
betaMin = 10; betaMax = 100;
mMin = 0.037; mMax = 0.52;
K = mMax/betaMin/phiMin*10;

params = cell(nSet,5);

for kk = 1:nSet
    r =   rMin + (rMax - rMin)*rand(nH,1);
    phi = phiMin + (phiMax - phiMin)*rand(1,nV);
    beta = betaMin + (betaMax - betaMin)*rand(1,nV);
    m  =   mMin + (mMax - mMin)*rand(1,nV);
    h = m./beta./phi;
    [h,idx] = sort(h,'descend');
    r = sort(r,'descend');
    m = m(idx)';
    phi = phi(idx);
    beta = beta(idx);
    beta = repmat(beta,nH,1);
    phi = repmat(phi,nH,1);    
    x0 = ( 0.95 + 0.1*rand(nH+nV,1)).*...
        equilibrium(matrices(:,:,focalM),r,a,K,phi,beta,m);
    params(kk,:) = {r,phi,beta,m,x0};  
end

clear kk r phi beta m h idx

save(strcat('data/',name))

%% parameters for the rule-based framework and low-nestedness net (fm1)
clear all
load('matrices/matrices_10by10_100_special_inv');

focalM = 1;
name = ['params_rb_fm' int2str(focalM) '_v2'];

nSet = 100;
[nH,nV] = size(matrices(:,:,1));
a = ones(nH);

rMin = 3.61; rMax=5/4*rMin; rMid = (rMin +rMax)/2;
mMax = 0.52; mMin = mMax*4/5; mMid = (mMin + mMax)/2;
phiMin = 2.4*10^-7; phiMax = 2.4*10^-6;
phiMid = (phiMin + phiMax)/2;
phi = phiMid*ones(1,nV);
beta = 30*ones(1,nV);
beta = repmat(beta,nH,1);
phi = repmat(phi,nH,1);
K = (2*mMid)/phiMid/30*10;% (m4+m7)/phi/beta<k

params = cell(nSet,5);

for kk = 1:nSet
    r = zeros(nH,1);
    r([ 1 3 5 9 2 6 8 10]) = rMid+ (rMax -rMid)*rand(8,1);
    r([4 7]) = rMin + (rMid - rMin)*rand(2,1);
    m = zeros(nH,1);
    m([ 1 3 5 9 2 6 8 10]) = mMid + (mMax - mMid)*rand(8,1);
    m([4 7]) = mMin + (mMid - mMin)*rand(2,1);
    x0 = ( 0.95 + 0.1*rand(nH+nV,1)).*...
    equilibrium(matrices(:,:,focalM),r,a,K,phi,beta,m);
    params(kk,:) = {r,phi,beta,m,x0}; 
end

clear kk r phi beta m h idx

save(strcat('data/',name))

%%
clear all
load('matrices/matrices_10by10_100_invertible');
focalM = 100;
deltaV = [0.005 0.01 0.02 0.05 0.1 0.2];
[nH,nV] = size(matrices(:,:,100));
nSet = 100;

a = ones(nH);
phiMin = 2.4*10^-8; phiMax = 2.4*10^-7;
betaMin = 10; betaMax = 100;
Hmin = 10^3; Hmax = 10^4;
Vmin = 10^6; Vmax = 10^7;
K = Hmax*10;

M = matrices(:,:,focalM);
phiStar = phiMin + (phiMax - phiMin)*rand(1,nV);
betaStar = betaMin + (betaMax - betaMin)*rand(1,nV);
betaStarM = repmat(betaStar,nH,1);
phiStarM = repmat(phiStar,nH,1);
Hstar = Hmin + (Hmax-Hmin)*rand(nH,1);
Vstar = Vmin + (Vmax - Vmin)*rand(nV,1);
mStar = (betaStarM.*phiStarM.*M)'*Hstar;
rStar = (phiStarM.*M)*Vstar./(1 - a*Hstar/K);


i = 1;

for delta = deltaV

    name = ['params_fm' num2str(focalM) '_delta' num2str(i) '_v4'];
    i = i+1;
    
    params = cell(nSet,5);

    for kk = 1:nSet
        r = rStar.*(1 - delta + 2*delta*rand(nH,1));
        m =  mStar.*(1 - delta + 2*delta*rand(nH,1));
        phi = phiStar.*(1 - delta + 2*delta*rand(1,nV)); 
        beta = betaStar.*(1 - delta + 2*delta*rand(1,nV));
        beta = repmat(beta,nH,1);
        phi = repmat(phi,nH,1);    
        x0 = ( 0.95 + 0.1*rand(nH+nV,1)).*[Hstar;Vstar];
        params(kk,:) = {r,phi,beta,m,x0};       
    end 

    clear kk r phi beta m h idx
    
    save(strcat('data/',name))
end

%%
clear all
load('matrices/matrices_10by10_100_invertible');
focalM = 1;
deltaV = [0.005 0.01 0.02 0.05 0.1 0.2];
[nH,nV] = size(matrices(:,:,100));
nSet = 100;

a = ones(nH);
phiMin = 2.4*10^-8; phiMax = 2.4*10^-7;
betaMin = 10; betaMax = 100;
Hmin = 10^3; Hmax = 10^4;
Vmin = 10^6; Vmax = 10^7;
K = Hmax*10;

M = matrices(:,:,focalM);
phiStar = phiMin + (phiMax - phiMin)*rand(1,nV);
betaStar = betaMin + (betaMax - betaMin)*rand(1,nV);
betaStarM = repmat(betaStar,nH,1);
phiStarM = repmat(phiStar,nH,1);
Hstar = Hmin + (Hmax-Hmin)*rand(nH,1);
Vstar = Vmin + (Vmax - Vmin)*rand(nV,1);
mStar = (betaStarM.*phiStarM.*M)'*Hstar;
rStar = (phiStarM.*M)*Vstar./(1 - a*Hstar/K);

i = 1;

for delta = deltaV

    name = ['params_fm' num2str(focalM) '_delta' num2str(i) '_v4'];
    i = i+1;
    
    params = cell(nSet,5);

    for kk = 1:nSet
        r = rStar.*(1 - delta + 2*delta*rand(nH,1));
        m =  mStar.*(1 - delta + 2*delta*rand(nH,1));
        phi = phiStar.*(1 - delta + 2*delta*rand(1,nV)); 
        beta = betaStar.*(1 - delta + 2*delta*rand(1,nV));
        beta = repmat(beta,nH,1);
        phi = repmat(phi,nH,1);    
        x0 = ( 0.95 + 0.1*rand(nH+nV,1)).*[Hstar;Vstar];
        params(kk,:) = {r,phi,beta,m,x0};       
    end 

    clear kk r phi beta m h idx
    
    save(strcat('data/',name))
end

%%
clear all
load('matrices/matrices_10by10_100_invertible');
focalM = 28;
deltaV = [0.005 0.01 0.02 0.05 0.1 0.2];
[nH,nV] = size(matrices(:,:,100));
nSet = 100;

a = ones(nH);
phiMin = 2.4*10^-8; phiMax = 2.4*10^-7;
betaMin = 10; betaMax = 100;
Hmin = 10^3; Hmax = 10^4;
Vmin = 10^6; Vmax = 10^7;
K = Hmax*10;

M = matrices(:,:,focalM);
phiStar = phiMin + (phiMax - phiMin)*rand(1,nV);
betaStar = betaMin + (betaMax - betaMin)*rand(1,nV);
betaStarM = repmat(betaStar,nH,1);
phiStarM = repmat(phiStar,nH,1);
Hstar = Hmin + (Hmax-Hmin)*rand(nH,1);
Vstar = Vmin + (Vmax - Vmin)*rand(nV,1);
mStar = (betaStarM.*phiStarM.*M)'*Hstar;
rStar = (phiStarM.*M)*Vstar./(1 - a*Hstar/K);

i = 1;

for delta = deltaV

    name = ['params_fm' num2str(focalM) '_delta' num2str(i) '_v4'];
    i = i+1;
    
    params = cell(nSet,5);

    for kk = 1:nSet
        r = rStar.*(1 - delta + 2*delta*rand(nH,1));
        m =  mStar.*(1 - delta + 2*delta*rand(nH,1));
        phi = phiStar.*(1 - delta + 2*delta*rand(1,nV)); 
        beta = betaStar.*(1 - delta + 2*delta*rand(1,nV));
        beta = repmat(beta,nH,1);
        phi = repmat(phi,nH,1);    
        x0 = ( 0.95 + 0.1*rand(nH+nV,1)).*[Hstar;Vstar];
        params(kk,:) = {r,phi,beta,m,x0};       
    end 

    clear kk r phi beta m h idx
    
    save(strcat('data/',name))
end