% Batch mode is run to test on different maxDrift, sampleName, and noise at 
% the same time. 
% To test on various maxDrift values, define the set "maxDriftAll" whose
% maximum is maxMaxDriftAll

%% Generate Model: 

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
[S0,L0,WGT,Lnormalizer] = generateModel(sampleName,[],WSz,NTheta,NTau,driftGT0,noiseLevel,0);

% using the same test image WGT, generate forward operator and 
% sinogram WITH drift error
rng('default');
driftGT1 = (2*rand(NTau, 1)- 1) * maxDrift;
driftGT1(1:ceil(maxMaxDriftAll)) = 0; 
driftGT1(end-ceil(maxMaxDriftAll)+1:end) = 0;
[S1,L1] = generateModel('',WGT,WSz,NTheta,NTau,driftGT1,noiseLevel,LNormalizer);


LAll = cell(2,1); SAll = LAll;
LAll{1} = L0; LAll{2} = L1;
SAll{1} = S0; SAll{2} = S1;

%% Show W and S
figNo = 10;
figure(figNo+1); imagesc(WGT); axis image; axis off; 
title('Object Ground Truth');

figure(figNo+2); imagesc(SAll{1}); axis image; axis off; 
title('Sinogram Without drift');

figure(figNo+3); imagesc(SAll{2}); axis image; axis off; 
title(sprintf('Sinogram With drift, PSNR=%.2fdB', difference(SAll{2}, SAll{1})));
tilefigs;

%% How good is calclating drift given the true sinogram pairs ??

interpMethods = {'linear', 'quadratic'}; 
driftComputed = calcDrift(SAll{2}, SAll{1}, maxDrift, interpMethods{1});
% zero out the ground truth for the begin and end empty beams. 
iStart = find(cumsum(driftComputed), 1, 'first');
iEnd = NTau + 1 - find(cumsum(driftComputed(end:-1:1)), 1, 'first');
idxHasDrift = iStart:iEnd;
idxNoDrift = [1:iStart-1, iEnd+1:NTau];

driftGT2 = driftGT1;         driftGT2(idxNoDrift) = 0;
figure(1001), clf; plot([driftGT2, driftComputed], '-x'); title(sprintf('psnr=%.2fdB', difference(driftComputed, driftGT2))); legend('GT', 'rec-testInterp'); grid on; grid minor; %driftComputed-driftGT2, 'dif'


%% acutual reconstruction with bad starting point (need to run recTompographyWithDrift first)
%drift = driftAll(:,end);    drift(idxNoDrift) = 0;
%figure(1002), clf; plot([driftGT2, drift, '-x']); title(sprintf('psnr=%.2fdB', difference(drift, driftGT2))); legend('GT', 'rec-actual'); %driftComputed-driftGT2, 'dif'

%%

%% II: Solve linear inverse problem: L*W = S
% solve for three cases. 
cases = {'with no drift', 'with known drift', 'with unknown drift'};
idxS = [1, 2, 2]; % corresponding S matrix for three cases
idxL = [1, 2, 1]; % corresponding L matrix for three cases, note in the last case, the drift is unkown, we use model of zero drift

%% (1) Reconstruction with solver from Wendy, just least square, no regularizer, not useful
%{
WWendy = zeros(Ny, Nx, numel(cases), 'like', WGT);
for n = 1:numel(cases)
    WWendy(:,:,n) = solveWendy(SAll(:,:,idxS(n)), LAll{idxL(n)}, WGT);
end
figNo = 30;
for n = 1:numel(cases)    
    figure(figNo+n); W = WWendy(:,:,n);  imagesc(W); axis image; title(sprintf('Rec Wendy %s, PSNR=%.2fdB', cases{n}, difference(W, WGT)));
end
tilefigs;
%}

%% Wendy-L2 regularizer. This is not very useful
%{
n = 3;
WWendyL2 = solveWendyL2Reg(SAll(:,:,idxS(n)), LAll{idxL(n)}, WGT);
figure(figNo+numel(cases)+1); W = WWendyL2;  imagesc(W); axis image; title(sprintf('Rec Wendy L2-regularizer %s, PSNR=%.2fdB', cases{n}, difference(W, WGT)));
%}

%% Reconstruction with solver from XH, with L1/TV regularizer.
% Need 100/500/1000+ iteration to get good/very good/excellent result with small regularizer.
% choose small maxSVDofA to make sure initial step size is not too small. 1.8e-6 and 1e-6 could make big difference for n=2 case
regType = 'TV'; % 'TV' or 'L1' % TV is better and cleaner for phantom example
maxIterA = 500; % 100 is not enough
regWt = 1e-5; % 1e-5 or 1e-8*max(WGT(:))/LNormalizer^2;  1e-6 to 1e-8 both good for phantom, use 1e-8 for brain, especically WGT is scaled to maximum of 1 not 40
maxSVDofA = 1e-4; % for TV, maxSVDofA = 1e-5/1e-4 or 1e-6/LNormalizer^2 NOT maxSVDofA = 1e-6/LNormalizer,  make sure it is actually much smaller than svds(L, 1), so that first step in TwIST is long enough 
paraTwist = {'xSz', WSz, 'regFun', regType, 'regWt', regWt, 'isNonNegative', 1, 'maxIterA', maxIterA, 'xGT', WGT, 'maxSVDofA', maxSVDofA, 'tolA', 1e-8};
WUncalib = zeros(WSz(1), WSz(2), numel(cases), 'like', WGT);

for n = 1:numel(cases)
    WUncalib(:,:,n) = solveTwist(SAll{idxS(n)}, LAll{idxL(n)}, paraTwist{:});
end
figNo = 50;
for n = 1:numel(cases)    
    figure(figNo+n); W = WUncalib(:,:,n);  imagesc(W); axis image; title(sprintf('Rec  %s, PSNR=%.2fdB, %s, regWt=%.1e, maxIter=%d', cases{n}, difference(W, WGT), regType, regWt, maxIterA));
end
tilefigs;

%% Try different regWt 
regWtAllUncalib = 1e-5*10.^(-3:3);
NWt = numel(regWtAllUncalib);
WUncalibAllWt = zeros(WSz(1), WSz(2), numel(cases), NWt, 'like', WGT);
psnrAllWtUncalib = zeros(numel(cases), NWt);
WUncalibBestWt = zeros(WSz(1), WSz(2), numel(cases), 'like', WGT);
regWtBestUncalib = zeros(numel(cases), 1);
for n = 1:numel(cases)
    for nWt = 1:NWt        
        paraTwistTry = {'xSz', WSz, 'regFun', regType, 'regWt', regWtAllUncalib(nWt), 'isNonNegative', 1, 'maxIterA', maxIterA, 'xGT', WGT, 'maxSVDofA', maxSVDofA, 'tolA', 1e-8};
        WUncalibAllWt(:,:,n,nWt) = solveTwist(SAll{idxS(n)}, LAll{idxL(n)}, paraTwistTry{:});
        psnrAllWtUncalib(n, nWt) = psnr(WUncalibAllWt(:,:,n,nWt), WGT, 1);  % WGT is normalized with maximum as 1
    end
    [~, idx] = max(psnrAllWtUncalib(n,:));
    regWtBestUncalib(n) = regWtAllUncalib(idx);
    WUncalibBestWt(:,:,n) = WUncalibAllWt(:,:,n,idx);    
end
%regWtBest = max(psnrAll
figNo = 100;
for n = 1:numel(cases)    
    figure(figNo+n); W = WUncalibBestWt(:,:,n);  imagesc(W); axis image; title(sprintf('Rec twist %s, PSNR=%.2fdB, %s, regWt=%.1e, maxIter=%d', cases{n}, difference(W, WGT), regType, regWtBestUncalib(n), maxIterA));
end
tilefigs;

% setup the regWt for our calibrated algorithm
regWt = regWtBestUncalib(3); % use regWtBest(2) where drift is known, or regWtBest(1) where drift is 0
paraTwist = {'xSz', WSz, 'regFun', regType, 'regWt', regWt, 'isNonNegative', 1, 'maxIterA', maxIterA, 'xGT', WGT, 'maxSVDofA', maxSVDofA, 'tolA', 1e-8};

%% Developping Stage: Two steps reconstruction, first TV on sinogram, then TV on object, just 1dB improvement for distorted one, not very useful
%Reconstruction with solver from XH, with L1/TV regularizer.
%{
regType = 'TV'; % 'TV' or 'L1' % TV is better and cleaner for phantom example
regWt = 1e-8*max(WGT(:)); % 1e-6 to 1e-8 both good for phantom, use 1e-8 for brain, especically WGT is scaled to maximum of 1 not 40
maxIterA = 500; % 100 is not enough
maxSVDofA = 1e-6; %svds(L, 1)*1e-4; % multiply by 1e-4 to make sure it is small enough so that first step in TwIST is long enough 
paraTwist = {'xSz', WSz, 'regFun', regType, 'regWt', regWt, 'isNonNegative', 1, 'maxIterA', maxIterA, 'xGT', WGT, 'maxSVDofA', maxSVDofA};

regWtSino = regWt*10000;
paraTwistSino = {'xSz', SSz, 'regFun', regType, 'regWt', regWtSino, 'isNonNegative', 1, 'maxIterA', maxIterA, 'xGT', SAll(:,:,1), 'maxSVDofA', maxSVDofA};

WUncalib = zeros(WSz(1), WSz(2), numel(cases), 'like', WGT);
SAllReg = zeros(size(SAll), 'like', SAll);


IdentityMat = speye(NTheta*NTau, NTheta*NTau);
for n = 1:size(SAll,3)
    SAllReg(:,:,n) = solveTwist(SAll(:,:,n), IdentityMat, paraTwistSino{:});  % may need to change TV just along one direction
end
for n = 1:numel(cases)
%for n = 1:1    
    WUncalib(:,:,n) = solveTwist(SAllReg(:,:,idxS(n)), LAll{idxL(n)}, paraTwist{:});
end
figNo = 60;
for n = 1:numel(cases)    
    figure(figNo+n); W = WUncalib(:,:,n);  imagesc(W); axis image; title(sprintf('Rec twist %s, PSNR=%.2fdB, %s, regWt=%.1e, maxIter=%d', cases{n}, difference(W, WGT), regType, regWt, maxIterA));
end
figNo = figNo+numel(cases);
sinoCases = {'no drift', 'drift'};
for n = 1:size(SAllReg,3)
    figure(figNo+n); SReg = SAllReg(:,:,n);  imagesc(SReg); axis image; title(sprintf('Sino %s smoothed, vs Sino No Drift PSNR=%.2fdB, %s, regWt=%.1e, maxIter=%d', sinoCases{n}, difference(SReg, SAll(:,:,1)), regType, regWtSino, maxIterA));
end

tilefigs([], [], [], 3); %tilefigs(handles,resize,nRows,nCols,...)
%}

%% Find the drift and Reconstruct by Alternative Minimization Method
% hack, just use script, and pass maxDrift here, to change to function
recTompographyWithDrift

%%
if ~batchMode && false
%driftIter = driftOld;
save('ret_20181121.mat', 'X', 'driftAll', 'driftGT', 'driftIter', 'maxDrift', 'SAll');
%%
paper_figures

%% How good is our interpolcation model, 
% even given exact drift, our interpolation model is not accurate enough to produce very good reconstruction. PSNR = 24/25dB 
% while use the exact model, PSNR = 34dB/50dB+ use Wendy or TwIST
driftGT = driftGTAll{2};
PGT =  genDriftMatrix(driftGT, NTheta, NTau);
L =  PGT*LAll{1};
%difference(L, LAll{2})  % psnr/psnrSparse = 30.16dB/9.63dB
difference(L*WGT(:), LAll{2}*WGT(:)) % psnr/psnrSparse = 34.82dB/32.31dB for phantom, but 45.34dB/44.80dB for brain which is more smooth
figNo = 70;
figure(figNo+1); W2 = solveTwist(SAll{2}, L, paraTwist{:});figure(figNo+4); W = W2;  imagesc(W); axis image; title(sprintf('Rec twist %s, PSNR=%.2fdB, %s, regWt=%.1e, maxIter=%d', 'Use Interpolated Forward Model', difference(W, WGT), regType, regWt, maxIterA));
figure(figNo+2); W2 = solveWendy(SAll{2}, L, WGT); figure(figNo+5); W = W2;  imagesc(W); axis image; title(sprintf('Rec Wendy %s, PSNR=%.2fdB', 'Use Interpolated Forward Model', difference(W, WGT)));

end