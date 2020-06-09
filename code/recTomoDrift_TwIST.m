function [W_history,S_history,info] = recTomoDrift_TwIST(WGT,L0,L1,S,maxIter,eta,...
                            maxDrift,epsilon,lambSet,driftGT,params)
                        
% convert the script recTomographyWithDrift.m to a function

if ~exist('maxDrift', 'var') || isempty(maxDrift)
    maxDrift = 1; 
end

% regType = params{4};
% regWt = params{6};
% maxIterA = params{10};



[Ny,Nx] = size(WGT);
[NTheta, NTau] = size(S);
WSz = size(WGT);


driftAll = zeros(NTau, maxIter); 
L = L0; 

WPrev = zeros(Ny, Nx); W = WPrev;
paraTwist2 = params;

scaleBegin = 1e2; %1e2 in the beginning(k=1), 
scaleEnd = 1e0; %1e0~1e1 in the end (k=maxIter),

psnr_history = zeros(maxIter,1);
W_history = zeros(Ny*Nx,maxIter);
S_history = zeros(NTheta*NTau,maxIter);
lambda_history = zeros(maxIter,1);
lambda_history(1) = paraTwist2{6};

regFun_history = zeros(maxIter,1);
objFun_history = zeros(maxIter,1);

etaeps_history = eta*0.01*norm(S(:))*ones(maxIter,1);    

% epsilon = noiseLevel*norm(S(:));


for k = 1:maxIter 

    if strcmp(lambSet,'old') == 1
        if maxIter > 1
            paraTwist2{6} = (1/scaleBegin)*params{6}*(scaleBegin + (k-1)*(scaleEnd-scaleBegin)/(maxIter-1)); 
        else
            paraTwist2{6} = (1/scaleBegin)*params{6}*scaleEnd;
        end
    end

%     etaeps_history(k) = etaeps_history(k)/(1.5^(k-1));

    eta = IRnormest(L1-L,S(:));
%     etaeps_history(k) = eta + epsilon;
    etaeps_history(k) = norm(epsilon + (L-L1)*WGT(:));
    lambda_history(k) = paraTwist2{6};
    
    [W,~,objFunVals,~,~,~,~,regFunVals] = solveTwist(S(:), L, 'x0', WPrev, paraTwist2{:}); 
    WPrev = W;
    
    regFun_history(k) = paraTwist2{6}*regFunVals(end);
    objFun_history(k) = objFunVals(end);
    alpha = objFunVals(end);
    beta = regFunVals(end);
    paraTwist2{6} = abs((etaeps_history(k) - alpha)/beta);
    

        
    S0 = reshape(L0*W(:), NTheta, NTau);
    drift =  calcDrift(S, S0, maxDrift);
    driftAll(:, k) = drift;
    
    psnr_history(k) = psnr(W, WGT);
    W_history(:,k) = W(:);
    S_history(:,k) = L*W(:);
            

    if k < maxIter
        P =  genDriftMatrix(drift, NTheta, NTau); L =  P*L0; % update L^k
        % L = XTM_Tensor_XH(WSz, NTheta, NTau, drift)/LNormalizer;
    end
    
    
    if k>1
        idxNoZero = 25:125;
        fprintf('iter=%d, mean squared error difference of new and old drift: %.2e\n', k, immse(driftAll(idxNoZero, k-1), driftAll(idxNoZero, k)));
        fprintf('iter=%d, mean squared error difference of new and GT drift: %.2e\n', k, immse(driftAll(idxNoZero, k), driftGT(idxNoZero, :)));
    end   
   
%     if k > 1 && psnr_history(k) < psnr_history(k-1)
%         psnr_history(k:end) = [];
%         W_history(:,k:end) = [];
%         S_history(:,k:end) = [];
%         break
%     end
        
    
end


% remove all the unused driftAll
% if k < size(driftAll, 2)
%     driftAll(:, k) = [];
% end



% show final result
% figNo = figNo + 1;
% figure(figNo);clf; figNo=figNo+1; plot([driftGT, drift]); legend('drift: GT', 'drift: Rec')
% SRec = reshape(L*W(:), NTheta, NTau); %%SRec2 = interpSino(S0, drift, 'linear'); difference(SRec2, SRec) % 317.24dB matches well
% figure(figNo);figNo=figNo+1; multAxes(@imagesc, {WGT, W}); multAxes(@title, {'Ground Truth', 'Reconstruction Final'}); linkAxesXYZLimColorView(); 
% figure(figNo);figNo=figNo+1; multAxes(@imagesc, {S, SRec}); multAxes(@title, {'Sinogram-Measured', sprintf('Sina-from rec, psnr=%.2fdB', difference(SRec,S))}); linkAxesXYZLimColorView(); 
% difference(S, SRec)


% Final reconstruction use L1-norm
% P =  genDriftMatrix(drift, NTheta, NTau); L =  P*L0;


% W = solveTwist(S, L, params{:}); 
% figure; imagesc(W); axis image; axis off
% title(sprintf('Rec twist %s, PSNR=%.2fdB, %s, regWt=%.1e, maxIter=%d', 'Use Interpolated Forward Model', difference(W, WGT), regType, regWt, maxIterA));
% W_history(:,end) = W(:);
% S_history(:,end) = L*W(:);
% 
% LExactfromDrift = XTM_Tensor_XH(WSz, NTheta, NTau, drift)/LNormalizer;
% 
% info.psnr_history = psnr_history;
% info.LExactfromDrift = LExactfromDrift;


info.psnr_history = psnr_history;
info.lambda_history = lambda_history;
info.regFun_history = regFun_history;
info.objFun_history = objFun_history;
info.etaeps_history = etaeps_history;
%etaeps