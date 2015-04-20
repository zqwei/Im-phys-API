%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 6.  Collected population decision decodability over time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison_V1_6

% addpath('../Func');
% setDir;
% 
% numFold             = 10;
% numTrials           = 200;
% numTestTrials       = 100;
% numTrainingTrials   = numTrials - numTestTrials;
% trainingTargets     = [true(numTrainingTrials/2,1); false(numTrainingTrials/2,1)];
% trainingTargets     = trainingTargets(randperm(numTrainingTrials));
% testTargets         = [true(numTestTrials/2,1); false(numTestTrials/2,1)];
% testTargets         = testTargets(randperm(numTestTrials));
% totTargets          = [testTargets; trainingTargets];
% 
% trainingDecisions   = trainingTargets(randperm(numTrainingTrials));
% testDecisions       = testTargets(randperm(numTestTrials));
% totDecisions        = [testDecisions; trainingDecisions];
% 
% load ([TempDatDir 'DataList.mat']);
% addNoise         = [1 1 0 0 0];
% 
% 
% for nData             = 4:4%length(DataSetList)
%     load([TempDatDir DataSetList(nData).name '.mat'])
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
%     xlabel('Time (s)');
%     ylabel('Decodability');
%     setPrint(4, 3, [PlotDir 'Collected_Units_Decodability__' DataSetList(nData).name], 'pdf')
% end

% for nData             = 4:4%length(DataSetList)
%     load([TempDatDir DataSetList(nData).name '.mat'])
%     figure;
%     nSessionData = shuffleSessionDataWithError(nDataSet, totTargets, totDecisions);
%     decodability = decodabilityLDA(nSessionData +randn(size(nSessionData))*1e-3/sqrt(numTrials)* addNoise(nData), trainingTargets, testTargets);
%     hold on
%     plot(DataSetList(nData).params.timeSeries, decodability, 'k', 'linewid',1);
%     xlim([min(DataSetList(nData).params.timeSeries) max(DataSetList(nData).params.timeSeries)]);
%     ylim([0.5 1])
%     gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 0.5)
%     box off;
%     hold off;
%     xlabel('Time (s)');
%     ylabel('Decodability');
%     setPrint(4, 3, [PlotDir 'Collected_Units_Decodability__' DataSetList(nData).name], 'pdf')
% end


addpath('../Func');
setDir;

numFold             = 30;
numTrials           = 200;
numTestTrials       = 100;
numTrainingTrials   = numTrials - numTestTrials;


load ([TempDatDir 'DataList.mat']);
addNoise         = [1 0 0 0 0];


for nData             = 1:length(DataSetList)
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
        
        nSessionData        = shuffleSessionData(nDataSet, totTargets, numTestTrials);
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
    setPrint(4, 3, [PlotDir 'Collected_Units_Decodability__' DataSetList(nData).name], 'pdf')
end
