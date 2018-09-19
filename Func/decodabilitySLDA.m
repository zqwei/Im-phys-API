%
% coeffSLDA.m
%
%
% ----------------------------
% Output:
%
% version 2.0
%
% based on paper:
%
% Guo et al., Biostatistics, 2007
%
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 


function decodability = decodabilitySLDA(nSessionData, trainingTargets, testTargets)
    
    T                 = size(nSessionData, 3);
    numTest           = length(testTargets);
    testData          = nSessionData(1:numTest, :, :);
    trainData         = nSessionData(numTest+1:end, :, :);

    decodability      = arrayfun(@(tIndex) 1 - loss(fitcdiscr...
                        (trainData(:,:,tIndex), trainingTargets, ...
                        'discrimType','pseudoLinear','Delta', 1, 'Gamma', 0.01),...
                        testData(:,:,tIndex), testTargets,...
                        'lossfun', 'classiferror'), ...
                        1:T, 'UniformOutput', false);
                    
    decodability      = cell2mat(decodability);
    
end