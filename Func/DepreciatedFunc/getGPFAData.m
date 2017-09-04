% 
% getGPFAData.m
% 
% version 1.0
%
% Comparison list
%
% Output:
% SpikeDataSet     --- yDim x nStim x nDecision (usually ignored) x T x N
%
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 

function dat     = getGPFAData(nDataSet, totTargets)
    nSessionData = shuffleSessionData(nDataSet, totTargets);
    lengthDat    = length(totTargets);
    dat          = repmat(struct('trialId',1, 'spikes', 1, ...
                                'trialType', 1), lengthDat, 1);
                            
    for nTrial = 1:lengthDat
        dat(nTrial).trialId                                             = nTrial;
        dat(nTrial).spikes                                              = squeeze(nSessionData(nTrial, :, :));
        dat(nTrial).trialType                                           = totTargets(nTrial);
    end

end