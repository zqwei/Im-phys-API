% 
% getGPFAData3D.m
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

function dat   = getGPFAData3D(nDataSet)

    lengthDat  = length(nDataSet.unit_yes_trial_index) + length(nDataSet.unit_no_trial_index);

    dat        = repmat(struct('trialId',1, 'spikes', 1, ...
                                'trialType', 1), lengthDat, 1);
                            
    for nTrial = 1:length(nDataSet.unit_yes_trial_index)
        dat(nTrial).trialId                                             = nTrial;
        dat(nTrial).spikes                                              = squeeze(nDataSet.unit_yes_trial(:, nTrial, :));
        dat(nTrial).trialType                                           = 1;
    end

    for nTrial = 1:length(nDataSet.unit_no_trial_index)
        dat(length(nDataSet.unit_yes_trial_index)+nTrial).trialId       = length(nDataSet.unit_yes_trial_index)+nTrial;
        dat(length(nDataSet.unit_yes_trial_index)+nTrial).spikes        = squeeze(nDataSet.unit_no_trial(:, nTrial, :));
        dat(length(nDataSet.unit_yes_trial_index)+nTrial).trialType     = 0;
    end

end