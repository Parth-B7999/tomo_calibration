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
    noiseLevel = 0;
end

% specify size of test image
Nx = 100; Ny = Nx; WSz = [Ny,Nx];

% define size of sinogram
NTheta = 45; 
NTau = ceil(sqrt(sum(WSz.^2))) + ceil(maxMaxDriftAll)*2; 
NTau = NTau + rem(NTau-Ny,2);

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


%% Show W and S
figNo = 10;
figure(figNo+1); imagesc(WGT); axis image; axis off; 
title('Object Ground Truth');

figure(figNo+2); imagesc(S0); axis image; axis off; 
title('Sinogram Without drift');

figure(figNo+3); imagesc(S1); axis image; axis off; 
title(sprintf('Sinogram With drift, PSNR=%.2fdB', difference(S1, S0)));
tilefigs;

%% How good is calclating drift given the true sinogram pairs ??

interpMethods = {'linear', 'quadratic'}; 
driftComputed = calcDrift(S1, S0, maxDrift, interpMethods{1});
% zero out the ground truth for the begin and end empty beams. 
iStart = find(cumsum(driftComputed), 1, 'first');
iEnd = NTau + 1 - find(cumsum(driftComputed(end:-1:1)), 1, 'first');
idxHasDrift = iStart:iEnd;
idxNoDrift = [1:iStart-1, iEnd+1:NTau];

driftGT2 = driftGT1;         driftGT2(idxNoDrift) = 0;
figure(1001), clf; plot([driftGT2, driftComputed], '-x'); title(sprintf('psnr=%.2fdB', difference(driftComputed, driftGT2))); legend('GT', 'rec-testInterp'); grid on; grid minor; %driftComputed-driftGT2, 'dif'


%% Reconstruction with solver from XH, with L1/TV regularizer.
% Need 100/500/1000+ iteration to get good/very good/excellent result with small regularizer.
% choose small maxSVDofA to make sure initial step size is not too small. 1.8e-6 and 1e-6 could make big difference for n=2 case
regType = 'TV'; % 'TV' or 'L1' % TV is better and cleaner for phantom example
maxIterA = 500; % 100 is not enough
regWt = 1e-5; % 1e-5 or 1e-8*max(WGT(:))/LNormalizer^2;  1e-6 to 1e-8 both good for phantom, use 1e-8 for brain, especically WGT is scaled to maximum of 1 not 40
maxSVDofA = 1e-4; % for TV, maxSVDofA = 1e-5/1e-4 or 1e-6/LNormalizer^2 NOT maxSVDofA = 1e-6/LNormalizer,  make sure it is actually much smaller than svds(L, 1), so that first step in TwIST is long enough 
paraTwist = {'xSz', WSz, 'regFun', regType, 'regWt', regWt, 'isNonNegative', 1, 'maxIterA', maxIterA, 'xGT', WGT, 'maxSVDofA', maxSVDofA, 'tolA', 1e-8};

WUncalib = solveTwist(S1, L0, paraTwist{:});

figure(20); imagesc(WUncalib); axis image; title(sprintf('Baseline, PSNR=%.2fdB, %s, regWt=%.1e, maxIter=%d', difference(WUncalib, WGT), regType, regWt, maxIterA));



%%
regWt = 1e-5; % use regWtBest(2) where drift is known, or regWtBest(1) where drift is 0
paraTwist = {'xSz', WSz, 'regFun', regType, 'regWt', regWt, 'isNonNegative', 1, 'maxIterA', maxIterA, 'xGT', WGT, 'maxSVDofA', maxSVDofA, 'tolA', 1e-8};

[W_history,S_history,infooo] = recTomoDrift_TwIST(WGT,L0,S1,20,...
                            maxDrift,LNormalizer,driftGT1,paraTwist);


%%
figNo = 500;
LExactfromDrift = infooo.LExactfromDrift;
W = solveTwist(S1, LExactfromDrift, paraTwist{:}); 
figure(figNo+115); imagesc(W); axis image; axis off
title(sprintf('Rec twist %s, PSNR=%.2fdB, %s, regWt=%.1e, maxIter=%d', 'Use Exact Forward Model from Drift', difference(W, WGT), regType, regWt, maxIterA));

paraTwist2 = paraTwist;
paraTwist2{6} = paraTwist2{6}*10; axis off
W = solveTwist(S1, LExactfromDrift, paraTwist2{:}); 
figure(figNo+116); imagesc(W); axis image; title(sprintf('Rec twist %s, PSNR=%.2fdB, %s, regWt=%.1e, maxIter=%d', 'Use Exact Forward Model from Drift', difference(W, WGT), regType, paraTwist2{6}, maxIterA));

