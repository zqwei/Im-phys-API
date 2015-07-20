%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 10. Similarity of PCA and LDA coefficient vectors as function of time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison_V1_10

addpath('../Func');
setDir;
load ([TempDatDir 'DataList.mat']);
addNoise         = [1 0 0 0 0 0];
numTrials           = 5000;
totTargets          = [true(numTrials/2,1); false(numTrials/2,1)];%rand(numTrials, 1) > 0.5;

if ~exist([PlotDir '/Collected_Units_PCA_LDA'],'dir')
    mkdir([PlotDir '/Collected_Units_PCA_LDA'])
end

for nData             = 1%fileToAnalysis
    load([TempDatDir DataSetList(nData).name '.mat'])
    figure;
    hold on
    nSessionData = shuffleSessionData(nDataSet, totTargets, numTrials);
    nSessionData = nSessionData + randn(size(nSessionData))*1e-3/sqrt(numTrials)* addNoise(nData);
    nSessionData = normalizationDim(nSessionData, 2);
    coeffs       = coeffLDA([nSessionData; nSessionData], totTargets, totTargets);
    imagesc(DataSetList(nData).params.timeSeries, DataSetList(nData).params.timeSeries, coeffs'*coeffs);
    xlim([min(DataSetList(nData).params.timeSeries) max(DataSetList(nData).params.timeSeries)]);
    ylim([min(DataSetList(nData).params.timeSeries) max(DataSetList(nData).params.timeSeries)]);
    caxis([-1 1]);
    axis xy;
    gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0], 'Color','k','Linestyle','--','linewid', 0.5);
    box off;
    hold off;
    xlabel('Time (s)')
    ylabel('Time (s)')
    setPrint(6, 4.5, [PlotDir 'Collected_Units_PCA_LDA/Similarity_LDA_coeffAll_' DataSetList(nData).name], 'tif')
end
close all
%%%


for nData             = fileToAnalysis
    load([TempDatDir DataSetList(nData).name '.mat'])
    figure;
    hold on
    nSessionData = shuffleSessionData(nDataSet, totTargets, numTrials);
    nSessionData = nSessionData + randn(size(nSessionData))*1e-3/sqrt(numTrials)* addNoise(nData);
    nSessionData = normalizationDim(nSessionData, 2);
    coeffs       = coeffPCA(nSessionData);
    imagesc(DataSetList(nData).params.timeSeries, DataSetList(nData).params.timeSeries, coeffs'*coeffs);
    xlim([min(DataSetList(nData).params.timeSeries) max(DataSetList(nData).params.timeSeries)]);
    ylim([min(DataSetList(nData).params.timeSeries) max(DataSetList(nData).params.timeSeries)]);
    caxis([0 1]);
    axis xy;
    gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0], 'Color','k','Linestyle','--','linewid', 0.5);
    box off;
    hold off;
    xlabel('Time (s)')
    ylabel('Time (s)')
    setPrint(6, 4.5, [PlotDir 'Collected_Units_PCA_LDA/Similarity_PCA_coeffAll_' DataSetList(nData).name], 'tif')
end
close all
    
%%%    



for nData             = fileToAnalysis
    load([TempDatDir DataSetList(nData).name '.mat'])
    figure;
    hold on
    nSessionData = shuffleSessionData(nDataSet, totTargets, numTrials);
    nSessionData = nSessionData + randn(size(nSessionData))*1e-3/sqrt(numTrials)* addNoise(nData);
    nSessionData = normalizationDim(nSessionData, 2);
    coeffPCAs    = coeffPCA(nSessionData);
    coeffLDAs    = coeffLDA([nSessionData; nSessionData], totTargets, totTargets);
    imagesc(DataSetList(nData).params.timeSeries, DataSetList(nData).params.timeSeries, abs(coeffPCAs'*coeffLDAs));
    xlim([min(DataSetList(nData).params.timeSeries) max(DataSetList(nData).params.timeSeries)]);
    ylim([min(DataSetList(nData).params.timeSeries) max(DataSetList(nData).params.timeSeries)]);
    caxis([0 1]);
    axis xy;
    gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0], 'Color','k','Linestyle','--','linewid', 0.5);
    box off;
    hold off;
    xlabel('Time (s)')
    ylabel('Time (s)')
    setPrint(6, 4.5, [PlotDir 'Collected_Units_PCA_LDA/Similarity_PCA_LDA_coeffAll_' DataSetList(nData).name], 'tif')
end

close all
