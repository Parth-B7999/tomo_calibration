% -------------------------------------------------------------------------
%  CONFIGURE YOUR MAIN SCRIPT NAME HERE
% -------------------------------------------------------------------------
mainScript = 'mainScript1.m';   % <-- change if your file is named differently
% -------------------------------------------------------------------------

sampleNames   = {'Phantom', 'Brain'};
gaussianSTDs  = [0.01 0.02];
maxDrifts     = [1    5];
lambda0s      = [1.3e-5 1e-8];           % values to drop into paraTwist{6}

for iS = 1:numel(sampleNames)
    for iN = 1:numel(gaussianSTDs)
        for iD = 1:numel(maxDrifts)
            for iL = 1:numel(lambda0s)

                %----------------------------------------------------------
                %  House-keeping
                %----------------------------------------------------------
                clearvars -except mainScript sampleNames gaussianSTDs ...
                                   maxDrifts lambda0s iS iN iD iL
                close all force

                %----------------------------------------------------------
                %  Variables consumed by the main script
                %----------------------------------------------------------
                batchMode   = true;                       
                sampleName  = sampleNames{iS};              
                gaussianSTD = gaussianSTDs(iN);           
                maxDrift    = maxDrifts(iD);               
                lambda0     = lambda0s(iL);                 

                fprintf('\n>>> RUN: %-7s | σ = %.3f | maxDrift = %d | λ₀ = %.1e\n', ...
                        sampleName, gaussianSTD, maxDrift, lambda0);

                %----------------------------------------------------------
                %  Call the main reconstruction script
                %----------------------------------------------------------
                run(mainScript);      

                %----------------------------------------------------------
                %  Optional: save whatever results you need
                %  Comment out or adapt as desired.
                %----------------------------------------------------------
                outFile = sprintf('result_%s_sigma%.3g_drift%d_lambda%.0e.mat', ...
                  sampleName, gaussianSTD, maxDrift, lambda0);

                % List every variable you might want …
                varsWanted = { ...
                    'WRecs','SRecs','drift','driftAll','info', ...            % NEW version
                    'WRecs_old','SRecs_old','drift_old','driftAll_old','info_old', ... % OLD version
                    'sampleName','gaussianSTD','maxDrift','lambda0' ...       % bookkeeping
                };
                
                % … keep only those that really exist (some are created only for
                % Type-I / Type-III runs, so we guard against “variable not found”)
                vars2save = varsWanted(cellfun(@(v) exist(v,'var')==1, varsWanted));
                
                if ~isempty(vars2save)
                    save(outFile, vars2save{:}, '-v7.3');
                else
                    warning('No result variables found to save for %s', outFile);
                end
            end
        end
    end
end