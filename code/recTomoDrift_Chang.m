function [W_history,S_history,info] = recTomoDrift_Chang(WGT,L0,S,maxIter,...
                            solver,maxDrift,LNormalizer,driftGT,params)
                        
% convert the script recTomographyWithDrift.m to a function

if ~exist('maxDrift', 'var') || isempty(maxDrift)
    maxDrift = 1; 
end


[Ny,Nx] = size(WGT);
[NTheta, NTau] = size(S);
WSz = size(WGT);


driftAll = zeros(NTau, maxIter); 
L = L0; 

WPrev = zeros(Ny*Nx,1);
W0 = WPrev;

scaleBegin = 1e2; %1e2 in the beginning(k=1), 
scaleEnd = 1e0; %1e0~1e1 in the end (k=maxIter),

psnr_history = zeros(maxIter,1);
W_history = zeros(Ny*Nx,maxIter+1);
S_history = zeros(NTheta*NTau,maxIter+1);

for k = 1:maxIter 
    fprintf('iteration = %d ', k)

%     if maxIter > 1
%         paraTwist2{6} = params{6}*(scaleBegin + (k-1)*(scaleEnd-scaleBegin)/(maxIter-1)); 
%     else
%         paraTwist2{6} = params{6}*scaleEnd;
%     end

    if strcmp(solver,'FLSQR-NNR')
        [WAll, EnrmAll] = flsqr_nnr(L, S(:), W0, WGT(:), params);
    elseif strcmp(solver,'IRN-LSQR-NNR')
        [WAll, ~,~, EnrmAll] = irn_lsqr_nnr(L, S(:), WGT(:), W0, optnnr);
    end
    
    [~,ind] = min(EnrmAll);
    W = WAll(:,ind);    
    WPrev = W;
    W = reshape(W,Ny,Nx);
    
    fprintf('best solution given by iteration %d of inner solver\n', ind)
        
    S0 = reshape(L0*W(:), NTheta, NTau);
    drift =  calcDrift(S, S0, maxDrift);
    driftAll(:, k) = drift;
    
    psnr_history(k) = psnr(W, WGT);
    W_history(:,k) = W(:);
    S_history(:,k) = L*W(:);

    
    P =  genDriftMatrix(drift, NTheta, NTau); L =  P*L0; % update L^k
    % L = XTM_Tensor_XH(WSz, NTheta, NTau, drift)/LNormalizer;
    
    
    
    if k>1
        idxNoZero = 25:125;
        fprintf('iter=%d, mean squared error difference of new and old drift: %.2e\n', k, immse(driftAll(idxNoZero, k-1), driftAll(idxNoZero, k)));
        fprintf('iter=%d, mean squared error difference of new and GT drift: %.2e\n', k, immse(driftAll(idxNoZero, k), driftGT(idxNoZero, :)));
    end    
    
end


% Final reconstruction use L1-norm
% P =  genDriftMatrix(drift, NTheta, NTau); L =  P*L0;
[WAll,EnrmAll] = flsqr_nnr(L0, S(:), WPrev, WGT(:), params);
[~,ind] = min(EnrmAll);
W = reshape(WAll(:,ind),Ny,Nx);
figure; imagesc(W); axis image; axis off
title(sprintf('%s PSNR=%.2fdB', 'Use Interpolated Forward Model', difference(W, WGT)));
W_history(:,end) = W(:);
S_history(:,end) = L*W(:);

LExactfromDrift = XTM_Tensor_XH(WSz, NTheta, NTau, drift)/LNormalizer;

info.psnr_history = psnr_history;
info.LExactfromDrift = LExactfromDrift;