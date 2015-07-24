%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 7.  Saturation of decodeability with number of neurons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison_V1_7

addpath('../Func');
setDir;
numTrials           = 400;
numTestTrials       = 200;
numTrainingTrials   = numTrials - numTestTrials;
trainingTargets     = [true(numTrainingTrials/2,1); false(numTrainingTrials/2,1)];
trainingTargets     = trainingTargets(randperm(numTrainingTrials));
testTargets         = [true(numTestTrials/2,1); false(numTestTrials/2,1)];
testTargets         = testTargets(randperm(numTestTrials));
totTargets          = [testTargets; trainingTargets];
load ([TempDatDir 'DataList.mat']);
addNoise            = [1 0 0 0 0 0];


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 7.1  Plot of LDA decorder
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for nData             = 1:length(DataSetList)
%     load([TempDatDir DataSetList(nData).name '.mat'])   
%     nSessionData = shuffleSessionData(nDataSet, totTargets);
%     numUnits     = length(nDataSet);
%     coeffs       = coeffLDA(nSessionData +randn(size(nSessionData))*1e-3/sqrt(numTrials)* addNoise(nData), trainingTargets, testTargets);
%     figure;
%     hold on
%     imagesc(DataSetList(nData).params.timeSeries, 1:numUnits, coeffs);
%     xlim([min(DataSetList(nData).params.timeSeries) max(DataSetList(nData).params.timeSeries)]);
%     ylim([1 numUnits])
%     colormap(french(128, 2));
%     caxis([-1/sqrt(numUnits) 1/sqrt(numUnits)]);
%     colorbar
%     axis xy;
%     gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 0.5);
%     box off;
%     hold off;
%     ylabel('Neuronal index')
%     xlabel('Time (s)')
%     setPrint(8, 6, [PlotDir 'Collected_Units_LDA_Coeffs_' DataSetList(nData).name], 'pdf')
% end
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 7.2  Plot of Sparsness using Kurtosis
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for nData             = 1:length(DataSetList)
%     load([TempDatDir DataSetList(nData).name '.mat'])  
%     nSessionData = shuffleSessionData(nDataSet, totTargets);
%     numUnits     = length(nDataSet);
%     coeffs       = coeffLDA(nSessionData +randn(size(nSessionData))*1e-3/sqrt(numTrials)* addNoise(nData), trainingTargets, testTargets);
%     figure;
% %     subplot(2,1,1)
% %     hold on
% %     plot(DataSetList(nData).params.timeSeries, mean(coeffs),'-b','linewid', 1.0);
% %     xlim([min(DataSetList(nData).params.timeSeries) max(DataSetList(nData).params.timeSeries)]);
% %     ylim([-0.1 0.1])
% %     gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 0.5);
% %     box off;
% %     hold off;
% %     
% %     subplot(2,1,2)
%     hold on
%     plot(DataSetList(nData).params.timeSeries, kurtosis(coeffs),'-k','linewid', 1.0);
%     xlim([min(DataSetList(nData).params.timeSeries) max(DataSetList(nData).params.timeSeries)]);
%     gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],3,  'Color','r','Linestyle','--','linewid', 0.5);
%     box off;
%     hold off;
%     ylabel('Kurtosis')
%     xlabel('Time (s)')
%     setPrint(8, 6, [PlotDir 'Collected_Units_LDA_Coeffs_Kurtosis_' DataSetList(nData).name], 'pdf')    
% end

if ~exist([PlotDir '/Collected_Units_LDA_Coeffs__KickOut'],'dir')
    mkdir([PlotDir '/Collected_Units_LDA_Coeffs__KickOut'])
end

perKickOut = 0:0.05:0.90;
numFold    = 30;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 7.3  Plot of Sparsness using Kicking-out neurons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for nData             = fileToAnalysis
    load([TempDatDir DataSetList(nData).name '.mat'])  
    nSessionData = shuffleSessionData(nDataSet, totTargets, numTestTrials);
    numUnits     = length(nDataSet);
    decodability = decodabilityLDAKickOut(nSessionData +randn(size(nSessionData))*1e-3/sqrt(numTrials)* addNoise(nData), trainingTargets, testTargets, perKickOut, numFold);    
    figure;
    hold on
    imagesc(DataSetList(nData).params.timeSeries, perKickOut*100, decodability.mean);
    xlim([min(DataSetList(nData).params.timeSeries) max(DataSetList(nData).params.timeSeries)]);
    ylim([0 90])
    caxis([0.5 1])
    colorbar
    axis xy;
    gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 0.5);
    box off;
    hold off;
    xlabel('Time (s)');
    ylabel('% KO neurons');
    setPrint(8, 6, [PlotDir 'Collected_Units_LDA_Coeffs__KickOut/Collected_Units_LDA_Coeffs__KickOut_' DataSetList(nData).name], 'pdf')
end

close all