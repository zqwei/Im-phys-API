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

function decodability = decodabilityLDA(nSessionData, trainingTargets, testTargets)
    
    T                 = size(nSessionData, 3);

    decodability      = arrayfun(@(tIndex) mean(testTargets == classify(...
                                           nSessionData(1:length(testTargets),:,tIndex), ...
                                           nSessionData(length(testTargets)+1:end,:,tIndex), ...
                                           trainingTargets)), 1:T, 'UniformOutput', false);
                                       
    decodability      = cell2mat(decodability);
    
end
