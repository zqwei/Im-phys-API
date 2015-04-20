%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 6.  Collected population decision decodability over time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison_V1_6


numTrials           = 3000;
numTestTrials       = 600;
numTrainingTrials   = numTrials - numTestTrials;
trainingTargets     = rand(numTrainingTrials, 1) > 0.5;
testTargets         = rand(numTestTrials, 1) > 0.5;
totTargets          = [testTargets; trainingTargets];

load ('TempDat/DataList.mat');
addNoise         = [1 1 0 0];

for nData             = 1:length(DataSetList)
    load(['TempDat/' DataSetList(nData).name '.mat'])    
    figure;
    nSessionData = shuffleSessionData(nDataSet, totTargets);
    decodability = decodabilityLDA(nSessionData +randn(size(nSessionData))*1e-3/sqrt(numTrials)* addNoise(nData), trainingTargets, testTargets);
    hold on
    plot(DataSetList(nData).params.timeSeries, decodability, 'k', 'linewid',1);
    xlim([min(DataSetList(nData).params.timeSeries) max(DataSetList(nData).params.timeSeries)]);
    ylim([0.5 1])
    gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 0.5)
    box off;
    hold off;
    setPrint(4, 3, ['Plot/Collected_Units_Decodability_' DataSetList(nData).name], 'png')
end