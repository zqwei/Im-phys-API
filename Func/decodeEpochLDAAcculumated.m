%
% decodeEpochLDAAcculumated.m
%
% based on decodeEpochLDA.m
% ----------------------------
% Output:
%
% version 1.0
%
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 

function decodability = decodeEpochLDAAcculumated(nSessionData, trainingEpoch, testEpoch, numT, numPeriods)
    

    groups            = classify(nSessionData(1:length(testEpoch),:), ...
                                nSessionData(length(testEpoch)+1:end,:), ...
                                trainingEpoch,'linear',ones(length(unique(trainingEpoch)),1)*1/(length(unique(trainingEpoch))));                                       
    groups            = reshape(groups, [], numT);
    
    decodability      = zeros(numPeriods, numT);
    
    for nPeriods      = 1:numPeriods
        decodability(nPeriods, :) = mean(groups == nPeriods);
    end
    
end
