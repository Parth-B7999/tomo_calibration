% Batch mode is run to test on different maxDrift, sampleName, and noise at 
% the same time. 
% To test on various maxDrift values, define the set "maxDriftAll" whose
% maximum is maxMaxDriftAll

%% Generate Model: 
clear

if ~exist('batchMode', 'var') || isempty(batchMode)
    batchMode = false; 
end

if ~batchMode 
    sampleName = 'Brain';
    maxDrift = 1; % maximum drift error relative to constant scan step size
    maxMaxDriftAll = maxDrift;
    noiseLevel = 0.01;
end

% specify size of test image
Nx = 100; Ny = Nx; WSz = [Ny,Nx];

% define size of sinogram
NTheta = 45; 
NTau = ceil(sqrt(sum(WSz.^2))) + ceil(maxMaxDriftAll)*2; 
NTau = NTau + rem(NTau-Ny,2);
% NTau = 144;

% generate test image, forward operator and sinogram WITHOUT drift error
driftGT0 = 0; 
[S0,L0,WGT,LNormalizer] = genModel(sampleName,[],WSz,NTheta,NTau,driftGT0,noiseLevel,0);

% using the same test image WGT, generate forward operator and 
% sinogram WITH drift error
rng('default');
driftGT1 = (2*rand(NTau, 1)- 1) * maxDrift;
driftGT1(1:ceil(maxMaxDriftAll)) = 0; 
driftGT1(end-ceil(maxMaxDriftAll)+1:end) = 0;
[S1,L1] = genModel('',WGT,WSz,NTheta,NTau,driftGT1,noiseLevel,LNormalizer);

%%
x0 = zeros(Nx*Ny,1);
optnnr.eta = 1.01;
optnnr.nl = noiseLevel;

%% FLSQR-NNR(v)
optnnr.p = 1;
optnnr.maxIt = 100;
optnnr.reg = 0;
optnnr.eta = 1.01;
optnnr.nl = noiseLevel;
optnnr.svdbasis = 'v';
optnnr.thr = 1e-3;
optnnr.gamma = 10^-10;
[X_fnnrv, Enrm_fnnrv] = flsqr_nnr(L0, S1(:), x0, WGT(:), optnnr);

% for i = 1:size(X_fnnrv,2)
%     psnr_fnnrv(i) = psnr(reshape(X_fnnrv(:,i),Ny,Nx),WGT);
% end
% [~,indPsnr] = max(psnr_fnnrv)
[~,ind] = min(Enrm_fnnrv)

Wfnnrv = X_fnnrv(:,ind);
Wfnnrv(Wfnnrv <0) = 0;
Wfnnrv = reshape(Wfnnrv,Ny,Nx);
figure,imagesc(Wfnnrv), axis image, axis off
title(sprintf('FLSQR-NNR(v), PSNR=%.2fdB, Iter %d', psnr(Wfnnrv,WGT), ind));



%% block coordinate descent
optnnr.reg = 0;
[Whis_fnnrv,Shis_fnnrv,psnr_fnnrv] = recTomoDrift_Chang(WGT,L0,S1,20,...
                            'FLSQR-NNR',maxDrift,LNormalizer,driftGT1,optnnr);

Wcalib_fnnrv = reshape(Whis_fnnrv(:,end),Ny,Nx);                       
figure,imagesc(Wcalib_fnnrv), axis image, axis off
title(sprintf('%s PSNR=%.2fdB', 'Use Interpolated Forward Model', psnr_fnnrv(end)));
%% FLSQR-NNR
optnnr.p = 1;
optnnr.maxIt = 100;
optnnr.reg = 0;
optnnr.eta = 1.01;
optnnr.nl = noiseLevel;
optnnr.svdbasis = 'x';
optnnr.thr = 1e-3;
optnnr.gamma = 10^-10;
x0 = zeros(Nx*Ny,1);
[X_fnnr, Enrm_fnnr] = flsqr_nnr(L0, S1(:), x0, WGT(:), optnnr);

% for i = 1:size(X_fnnr,2)
%     psnr_fnnr(i) = psnr(reshape(X_fnnr(:,i),Ny,Nx),WGT);
% end
% [~,indPsnr] = max(psnr_fnnr)

[~,ind] = min(Enrm_fnnr)

Wfnnr = X_fnnr(:,ind);
Wfnnr(Wfnnr <0) = 0;
Wfnnr = reshape(Wfnnr,Ny,Nx);
figure,imagesc(Wfnnr), axis image, axis off
title(sprintf('FLSQR-NNR, PSNR=%.2fdB, Iter %d', psnr(Wfnnr,WGT), ind));

%% IRN-LSQR-NNR
optnnr.p = 2;
optnnr.gamma = 10^-10;
optnnr.maxIt = 25;
optnnr.cycles = 4;
optnnr.thr = 1e-3;
optnnr.reg = 0;
optnnr.weigthtype = 'sqrt'; 
optnnr.thrstop = 1e-8;

[X_irn_tot, ~,~, Enrm_irn] = irn_lsqr_nnr(L0, S1(:), WGT(:), x0, optnnr);

% for i = 1:size(X_irn_tot,2)
%     psnr_irn(i) = psnr(reshape(X_irn_tot(:,i),Ny,Nx),WGT);
% end
% [~,indPsnr] = max(psnr_irn)
[~,ind] = min(Enrm_irn)

Wirn = X_irn_tot(:,ind);
Wirn(Wirn <0) = 0;
Wirn = reshape(Wirn,Ny,Nx);
figure,imagesc(Wirn), axis image, axis off
title(sprintf('IRN-LSQR-NNR, PSNR=%.2fdB, Iter %d', psnr(Wirn,WGT), ind));

%% 
[Whis_irn,Shis_irn,info_irn] = recTomoDrift_Chang(WGT,L0,S1,10,...
                            'IRN-LSQR-NNR',maxDrift,LNormalizer,driftGT1,optnnr);
