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
% kfoldLoss function gives the error rate, aka percentage of misclassied error
% -------------------------------------------------------------------------
%
% version 3.0
% An explict control of training and test datasets are applied
% LDA based on fitcdiscr function
% loss function gives the error rate, aka percentage of misclassied error
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 

%%% 
%%% version 1.0
%%% 
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

%%% 
%%% version 2.0
%%% 
% function decodability = decodabilityLDA(nSessionData, kFold)
%     
%     T                 = size(nSessionData, 3);
% 
%     decodability      = arrayfun(@(tIndex) 1 - kfoldLoss(fitcdiscr...
%                         (nSessionData(:,:,tIndex), totTargets, ...
%                         'discrimType','pseudoLinear', 'KFold', kFold),...
%                         'lossfun', 'classiferror'), ...
%                         1:T, 'UniformOutput', false);
%                     
%     decodability      = cell2mat(decodability);
%     
% end


%%% 
%%% version 3.0
%%% 
function decodability = decodabilityLDA(nSessionData, trainingTargets, testTargets)
    
    T                 = size(nSessionData, 3);
    numTest           = length(testTargets);
    testData          = nSessionData(1:numTest, :, :);
    trainData         = nSessionData(numTest+1:end, :, :);

    decodability      = arrayfun(@(tIndex) 1 - loss(fitcdiscr...
                        (trainData(:,:,tIndex), trainingTargets, ...
                        'discrimType','pseudoLinear'),...
                        testData(:,:,tIndex), testTargets,...
                        'lossfun', 'classiferror'), ...
                        1:T, 'UniformOutput', false);
                    
    decodability      = cell2mat(decodability);
    
end