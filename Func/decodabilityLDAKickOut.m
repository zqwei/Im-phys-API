%
% decodabilityLDAKickOut.m
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

function decodability = decodabilityLDAKickOut(nSessionData, trainingTargets, testTargets, perKickOut, numFold)
    
    T                 = size(nSessionData, 3);
    numUnit           = size(nSessionData, 2);
    numKick           = length(perKickOut);
    
    decodabilityAllFold = zeros(numKick, T, numFold);

    for nFold         = 1:numFold
        
        decodabilityFold = arrayfun(@(tPer) decodabilityLDA(...
                                    nSessionData(:,rand(numUnit,1)>tPer,:)...
                                    , trainingTargets, testTargets), perKickOut, 'UniformOutput', false);
        decodabilityFold = cell2mat(decodabilityFold');
        decodabilityAllFold(:, :, nFold) = decodabilityFold;        
    end
    
    decodability.mean    = mean(decodabilityAllFold, 3);
    decodability.std     = std(decodabilityAllFold, [], 3)/sqrt(numFold);
    
end
