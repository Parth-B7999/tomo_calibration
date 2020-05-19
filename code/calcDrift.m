function drift = calcDrift(S2, S, maxDrift, interpMethod)
%compute the drift error from sinogram matching
% Input: 
%   S2:  the measured Sinogram with beamline drifted, NTheta*NTau 2D matrix, assume it is interploation from S here. 
%   S: the estimated Sinogram with beamline NOT drifted, NTheta*NTau 2D matrix
%       We estimate it as reshape(L0*X(:), NTheta, NTau), where X is the reconstructed object
%   maxDrift: the maximum possible drift amount in the unit of beamlines, integer here (or pre-convert to integer: not tested)
%   interpMethod:    interplolation options
%               'linear'
%               'quadratic', 
% Output: 
%   drift: the drift amount for each beamline (assume invariant with rotation angle),  NTau*1 vector
% See also
%   interpSino(); genDriftMatrix();

if nargin < 1; testMe(); end


if ~exist('maxDrift', 'var') || isempty(maxDrift)
    maxDrift = 1; 
end

if ~exist('interpMethod', 'var') || isempty(interpMethod)
    interpMethod = 'linear';
end

assert(maxDrift==floor(maxDrift), 'maxDrift need to be integer!, but it is %.5f', maxDrift)


%%
[NTheta, NTau] = size(S2);
SDiffFwd = [S(:, 2:end) - S(:, 1:end-1), zeros(NTheta, 1)]; % SDiffFwd(:,n) = S(:,n+1) - S(:,n); Foward difference

% drift = d+alpha, where d is the integer part, and 0 <= alpha <1 is the fractional part
d = zeros(NTau, 1);
alpha = zeros(NTau, 1);

switch lower(interpMethod)
    case 'linear'
        for i = maxDrift+1 : NTau-maxDrift % compute d(i), note zero boundaries d(i)=0 for i < maxDrift+1 or i > NTau-maxDrift
            if norm(S2(:,i), 2) < NTheta*eps && norm(S(:,i), 2) < NTheta*eps; continue;  end  % empty beam matched to empty beam, skip & set drift to 0
            objOpt = inf;
            for k = -maxDrift : maxDrift-1  % interpolate between i+k and i+k+1
                % find the optimal alpha(i)
                a = SDiffFwd(:, i+k);
                b = S2(:,i) - S(:, i+k);
                aTa = a'*a;
                aTb = a'*b;
                bTb = b'*b;
                if aTa < NTheta*eps
                    alphaCur = 0;
                else
                    alphaCur = aTb / aTa;
                end
                alphaCur = min(max(alphaCur, 0), 1); % clip to [0, 1]
                objCur = aTa*alphaCur^2 - 2*aTb*alphaCur + bTb;
                if objCur <= objOpt
                    objOpt = objCur;
                    d(i) = k;
                    alpha(i) = alphaCur;
                end
            end
        end
    case 'quadratic'
        % don't loop through all k, but rather compute and store result for all posible k, and pickup the best
        printf('todo\n');
    otherwise
        error('Unknown interpolation method: %s.', interpMethod);         
end


% return the total drift
drift = d + alpha;


end

function testMe()
load('ret_20181121.mat')
drift = calcDrift(SAll(:,:,2), SAll(:,:,1), maxDrift, 'linear');
%%
%interpMethods = {'linear', 'quadratic', 'linearRaw', 'quadraticRaw', 'linearOld'}; % [34.82, 34.99, 30.41, 34.74, 34.82]dB
interpMethods = {'linear'}; % [34.82, 34.99, 30.41, 34.74, 34.82]dB
driftAll = zeros(numel(driftGT), numel(interpMethods));
for n = 1:numel(interpMethods)
    driftAll(:,n) = calcDrift(SAll(:,:,2), SAll(:,:,1), maxDrift, interpMethods{n});
end

%% Clean the left/right border, set to 0
drift = calcDrift(SAll(:,:,2), SAll(:,:,1), maxDrift, 'linear');
driftGT2 = driftGT;
iStart = find(cumsum(drift), 1, 'first');
iEnd = numel(driftGT) + 1 - find(cumsum(drift(end:-1:1)), 1, 'first');
driftGT2([1:iStart-1, iEnd+1:end]) = 0;
%%
for n = 1:numel(interpMethods)    
    fprintf('Interpolation method: %s, psnr = %.2fdB.\n', interpMethods{n}, difference(driftAll(:,n),  driftGT2));   
end
%%
figure, plot([driftGT2, drift], '-x'); title(sprintf('psnr=%.2fdB', difference(drift, driftGT2))); legend('GT', interpMethods{1}); %driftComputed-driftGT2, 'dif'
driftIter2 = driftIter; driftIter2([1:iStart-1, iEnd+1:end]) = 0;
figure, plot([driftGT2, driftIter2], '-x'); title(sprintf('psnr=%.2fdB', difference(driftIter2, driftGT2))); legend('GT', [interpMethods{1}, '\_iter']); %driftComputed-driftGT2, 'dif'

end
