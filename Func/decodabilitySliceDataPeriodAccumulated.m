%
% decodabilitySliceDataPeriodAccumulated.m
%
% based on decodabilitySliceDataPeriod.m
% ----------------------------
% Output:
%
% version 1.0
%
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 

function decodability               = decodabilitySliceDataPeriodAccumulated(nSessionData, EpochIndex, sliceIndex, noiseData, numTestTrials, numPeriods)
    [numTrials, ~, numUnits]        = size(nSessionData);
    nSessionData                    = nSessionData(:,sliceIndex,:);
    EpochIndex                      = EpochIndex(:,sliceIndex);
    numT                            = length(sliceIndex);
    nIndex                          = randperm(numTrials);
    nTestData                       = reshape(nSessionData(nIndex(1:numTestTrials),:,:), numTestTrials*numT, numUnits);
    nTrainData                      = reshape(nSessionData(nIndex(numTestTrials+1:end),:,:), (numTrials-numTestTrials)*numT, numUnits);
    nSessionData                    = [nTestData; nTrainData];
    testEpoch                       = reshape(EpochIndex(1:numTestTrials,:), numTestTrials*numT, 1);
    trainingEpoch                   = reshape(EpochIndex(numTestTrials+1:end,:), (numTrials-numTestTrials)*numT, 1);
    decodability                    = decodeEpochLDAAcculumated(nSessionData +randn(size(nSessionData))*1e-3/sqrt(numTrials)* noiseData, trainingEpoch, testEpoch, numT, numPeriods);    
end
