%
% decodabilitySliceDataTaper.m
%
%
% ----------------------------
% Output:
%
% version 1.0
%
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 


function decodability               = decodabilitySliceDataTaper(nSessionData, EpochIndex, taperLength, noiseData, numTestTrials)
    [numTrials, numT, numUnits]     = size(nSessionData);
    nSessionData                    = nSessionData(:,taperLength+1:end,:);
    EpochIndex                      = EpochIndex(:,1:end-taperLength);    
    numT                            = numT - taperLength;
    nTestData                       = reshape(nSessionData(1:numTestTrials,:,:), numTestTrials*numT, numUnits);
    nTrainData                      = reshape(nSessionData(numTestTrials+1:end,:,:), (numTrials-numTestTrials)*numT, numUnits);
    nSessionData                    = [nTestData; nTrainData];
    testEpoch                       = reshape(EpochIndex(1:numTestTrials,:), numTestTrials*numT, 1);
    trainingEpoch                   = reshape(EpochIndex(numTestTrials+1:end,:), (numTrials-numTestTrials)*numT, 1);
    decodability                    = decodeEpochLDA(nSessionData +randn(size(nSessionData))*1e-3/sqrt(numTrials)* noiseData, trainingEpoch, testEpoch, numT);    
end