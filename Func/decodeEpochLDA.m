%
% decodabilityLDA.m
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

function decodability = decodeEpochLDA(nSessionData, trainingEpoch, testEpoch, numT)
    
    decodability      = testEpoch == classify(nSessionData(1:length(testEpoch),:), ...
                                                nSessionData(length(testEpoch)+1:end,:), ...
                                                trainingEpoch,'linear',ones(length(unique(trainingEpoch)),1)*1/(length(unique(trainingEpoch))));
                                       
    decodability      = mean(reshape(decodability, [], numT));
end
