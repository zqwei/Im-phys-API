%
% decodabilityLDA.m
%
%
% ----------------------------
% Output:
%
% version 1.0
% LDA based on classify function
% -------------------------------------------------------------------------
%
% version 2.0
% LDA based on fitcdiscr function
% loss function gives the error rate, aka percentage of misclassied error
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 

% function decodability = decodabilityLDA(nSessionData, trainingTargets, testTargets)
%     
%     T                 = size(nSessionData, 3);
% 
%     decodability      = arrayfun(@(tIndex) mean(testTargets == classify(...
%                                            nSessionData(1:length(testTargets),:,tIndex), ...
%                                            nSessionData(length(testTargets)+1:end,:,tIndex), ...
%                                            trainingTargets)), 1:T, 'UniformOutput', false);
%                                        
%     decodability      = cell2mat(decodability);
%     
% end

function decodability = decodabilityLDA(nSessionData, trainingTargets, testTargets)
    
    T                 = size(nSessionData, 3);

    decodability      = arrayfun(@(tIndex) 1 - loss(fitcdiscr...
                        (nSessionData(:,:,tIndex), totTargets, ...
                        'discrimType','pseudoLinear', 'KFold', kFold),...
                        'lossfun', 'classiferror'), ...
                        1:T, 'UniformOutput', false);
                    
    decodability      = cell2mat(decodability);
    
end