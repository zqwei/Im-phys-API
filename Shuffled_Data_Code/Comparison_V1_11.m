%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 11. DPCA etc.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison_V1_11

function Comparison_V1_11

    numTrials           = 3000;
    numTestTrials       = 600;
    numTrainingTrials   = numTrials - numTestTrials;
    trainingTargets     = rand(numTrainingTrials, 1) > 0.5;
    testTargets         = rand(numTestTrials, 1) > 0.5;
    totTargets          = [testTargets; trainingTargets];
    load ('TempDat/DataList.mat');
    addNoise            = [1 1 0 0];        
    % Comparison_V1_10_1(DataSetList, numTestTrials, totTargets, addNoise, numTrials)
    % Comparison_V1_10_2(DataSetList, numTestTrials, totTargets, addNoise, numTrials)
    % Comparison_V1_10_3(DataSetList, numTestTrials, totTargets, addNoise, numTrials)
    % Comparison_V1_10_4(DataSetList, numTestTrials, totTargets, addNoise, numTrials)
end