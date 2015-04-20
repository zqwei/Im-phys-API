%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 7.  Saturation of decodeability with number of neurons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison_V1_7


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 7.1  Plot of LDA decorder
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    nSessionData = shuffleSessionData(nDataSet, totTargets);
    numUnits     = length(nDataSet);
    coeffs       = coeffLDA(nSessionData +randn(size(nSessionData))*1e-3/sqrt(numTrials)* addNoise(nData), trainingTargets, testTargets);
    figure;
    hold on
    imagesc(DataSetList(nData).params.timeSeries, 1:numUnits, coeffs);
    xlim([min(DataSetList(nData).params.timeSeries) max(DataSetList(nData).params.timeSeries)]);
    ylim([1 numUnits])
    colormap(french(128, 2));
    caxis([-1/sqrt(numUnits) 1/sqrt(numUnits)]);
    colorbar
    axis xy;
    gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[], 'Color','w','Linestyle','--','linewid', 0.5);
    box off;
    hold off;
    setPrint(4, 3, ['Plot/Collected_Units_LDA_Coeffs_' DataSetList(nData).name], 'png')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 7.2  Plot of Sparsness using Kurtosis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numTrials           = 3000;
numTestTrials       = 600;
numTrainingTrials   = numTrials - numTestTrials;
trainingTargets     = rand(numTrainingTrials, 1) > 0.5;
testTargets         = rand(numTestTrials, 1) > 0.5;
totTargets          = [testTargets; trainingTargets];

load ('TempDat/DataList.mat');
addNoise            = [1 1 0 0];

for nData             = 1:length(DataSetList)
    load(['TempDat/' DataSetList(nData).name '.mat'])    
    nSessionData = shuffleSessionData(nDataSet, totTargets);
    numUnits     = length(nDataSet);
    coeffs       = coeffLDA(nSessionData +randn(size(nSessionData))*1e-3/sqrt(numTrials)* addNoise(nData), trainingTargets, testTargets);
    figure;
    subplot(2,1,1)
    hold on
    plot(DataSetList(nData).params.timeSeries, mean(coeffs),'-b','linewid', 1.0);
    xlim([min(DataSetList(nData).params.timeSeries) max(DataSetList(nData).params.timeSeries)]);
    ylim([-0.1 0.1])
    gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 0.5);
    box off;
    hold off;
    
    subplot(2,1,2)
    hold on
    semilogy(DataSetList(nData).params.timeSeries, kurtosis(coeffs),'-b','linewid', 1.0);
    xlim([min(DataSetList(nData).params.timeSeries) max(DataSetList(nData).params.timeSeries)]);
    gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],3,  'Color','k','Linestyle','--','linewid', 0.5);
    box off;
    hold off;
    setPrint(8, 6, ['Plot/Collected_Units_LDA_Coeffs_Kurtosis_' DataSetList(nData).name], 'png')    
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 7.3  Plot of Sparsness using Kicking-out neurons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numTrials           = 3000;
numTestTrials       = 600;
numTrainingTrials   = numTrials - numTestTrials;
trainingTargets     = rand(numTrainingTrials, 1) > 0.5;
testTargets         = rand(numTestTrials, 1) > 0.5;
totTargets          = [testTargets; trainingTargets];
perKickOut          = 0:0.05:0.95;
numFold             = 10;
load ('TempDat/DataList.mat');
addNoise            = [1 1 0 0];

for nData             = 1:length(DataSetList)
    load(['TempDat/' DataSetList(nData).name '.mat'])    
    nSessionData = shuffleSessionData(nDataSet, totTargets);
    numUnits     = length(nDataSet);
    decodability = decodabilityLDAKickOut(nSessionData +randn(size(nSessionData))*1e-3/sqrt(numTrials)* addNoise(nData), trainingTargets, testTargets, perKickOut, numFold);    
    figure;
    hold on
    imagesc(DataSetList(nData).params.timeSeries, perKickOut*100, decodability.mean);
    xlim([min(DataSetList(nData).params.timeSeries) max(DataSetList(nData).params.timeSeries)]);
    ylim([0 95])
    caxis([0.5 1])
    colorbar
    axis xy;
    gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 0.5);
    box off;
    hold off;
    setPrint(4, 3, ['Plot/Collected_Units_LDA_Coeffs_KickOut_' DataSetList(nData).name], 'png')
end
