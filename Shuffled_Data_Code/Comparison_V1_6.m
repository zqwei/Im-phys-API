%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 6.  Collected population decision decodability over time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison_V1_6


% numTrials           = 3000;
% numTestTrials       = 600;
% numTrainingTrials   = numTrials - numTestTrials;
% trainingTargets     = rand(numTrainingTrials, 1) > 0.5;
% testTargets         = rand(numTestTrials, 1) > 0.5;
% totTargets          = [testTargets; trainingTargets];
% 
% load ('TempDat/DataList.mat');
% addNoise         = [1 1 0 0];
% 
% for nData             = 1:length(DataSetList)
%     load(['TempDat/' DataSetList(nData).name '.mat'])    
%     figure;
%     nSessionData = shuffleSessionData(nDataSet, totTargets);
%     decodability = decodabilityLDA(nSessionData +randn(size(nSessionData))*1e-3/sqrt(numTrials)* addNoise(nData), trainingTargets, testTargets);
%     hold on
%     plot(DataSetList(nData).params.timeSeries, decodability, 'k', 'linewid',1);
%     xlim([min(DataSetList(nData).params.timeSeries) max(DataSetList(nData).params.timeSeries)]);
%     ylim([0.5 1])
%     gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 0.5)
%     box off;
%     hold off;
%     setPrint(4, 3, ['Plot/Collected_Units_Decodability_' DataSetList(nData).name], 'png')
% end


addpath('../Func');
setDir;

numFold             = 50;
numTrials           = 4000;
numTestTrials       = 2000;
numTrainingTrials   = numTrials - numTestTrials;
% numRandPickUnits    = 60;

load ([TempDatDir 'DataListShuffle.mat']);
numRandPickUnits    = 1000;
addNoise            = [1 0 0 0 0 0];
fileToAnalysis      = [1 5];

if ~exist([PlotDir '/Collected_Units_Decodability'],'dir')
    mkdir([PlotDir '/Collected_Units_Decodability'])
end

for nData             = fileToAnalysis
    load([TempDatDir DataSetList(nData).name '.mat'])
    figure;
    
    decodability = zeros(numFold, size(nDataSet(1).unit_yes_trial,2));
    
    for nFold    = 1:numFold
        trainingTargets     = [true(numTrainingTrials/2,1); false(numTrainingTrials/2,1)];
        trainingTargets     = trainingTargets(randperm(numTrainingTrials));
        testTargets         = [true(numTestTrials/2,1); false(numTestTrials/2,1)];
        testTargets         = testTargets(randperm(numTestTrials));
        totTargets          = [testTargets; trainingTargets];

        trainingDecisions   = trainingTargets(randperm(numTrainingTrials));
        testDecisions       = testTargets(randperm(numTestTrials));
        totDecisions        = [testDecisions; trainingDecisions];
        
        randPickUnits       = randperm(length(nDataSet));
        randPickUnits       = randPickUnits(1:numRandPickUnits);
                
        nSessionData        = shuffleSessionData(nDataSet(randPickUnits), totTargets, numTestTrials);
        decodability(nFold,:) = decodabilityLDA(nSessionData +randn(size(nSessionData))*1e-3/sqrt(numTrials)* addNoise(nData), trainingTargets, testTargets);
    end
    hold on
    plot(DataSetList(nData).params.timeSeries, mean(decodability,1), 'k', 'linewid',1);
    xlim([min(DataSetList(nData).params.timeSeries) max(DataSetList(nData).params.timeSeries)]);
    ylim([0.5 1])
    gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 0.5)
    box off;
    hold off;
    xlabel('Time (s)');
    ylabel('Decodability');
    setPrint(4, 3, [PlotDir 'Collected_Units_Decodability/All_Units_Decodability_FixedNumberUnits_' DataSetList(nData).name], 'pdf')
end
close all;