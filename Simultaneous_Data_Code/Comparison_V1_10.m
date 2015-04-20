%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 10. Similarity of PCA and LDA coefficient vectors as function of time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison_V1_10

addpath('../Func');
setDir;
load ([TempDatDir 'DataList.mat']);
addNoise         = [1 1 0 0 0];
numTrials           = 3000;
numTestTrials       = 600;
numTrainingTrials   = numTrials - numTestTrials;
trainingTargets     = rand(numTrainingTrials, 1) > 0.5;
testTargets         = rand(numTestTrials, 1) > 0.5;
totTargets          = [testTargets; trainingTargets];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 10.1 Fraction of variance captured by PCA, per session, by trial period.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% for nData             = 1:length(DataSetList)
%     load([TempDatDir DataSetList(nData).name '.mat'])   
%     nSessionData = shuffleSessionData(nDataSet, totTargets);
%     nSessionData = normalizationDim(nSessionData, 2);
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
%     setPrint(8, 6, ['Collected_Units_LDA_Coeffs_' DataSetList(nData).name], 'tif')
% end

for nData             = 1:length(DataSetList)
    load([TempDatDir DataSetList(nData).name '.mat'])
    numSession       = length(nDataSet3D);
    m                = ceil((numSession)/4);
    numUnit          = arrayfun(@(x) length(x.nUnit), nDataSet3D, 'UniformOutput', false);
    [~, sortUnit]    = sort([numUnit{:}],'descend');    
    nDataSet3D       = nDataSet3D(sortUnit);
    figure;
    for nPlot        = 1:numSession
        subplot(m, 4, nPlot)
        hold on
        nSessionData = [permute(nDataSet3D(nPlot).unit_yes_trial,[2, 1, 3]); permute(nDataSet3D(nPlot).unit_no_trial,[2, 1, 3])];
        nSessionData = normalizationDim(nSessionData, 2);
        nTargets     = [true(size(nDataSet3D(nPlot).unit_yes_trial, 2),1); false(size(nDataSet3D(nPlot).unit_no_trial, 2),1)];
        nSessionData = [nSessionData; nSessionData]; %#ok<AGROW>
        coeffs       = coeffLDA(nSessionData +randn(size(nSessionData))*1e-3/sqrt(numTrials)* addNoise(nData), nTargets, nTargets);
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
        title([num2str(length(nDataSet3D(nPlot).nUnit)) ' Neurons'])
    end

%     subplot(m, 4, nPlot+1)
%     caxis([-1 1])
%     colorbar
%     axis off
    setPrint(4*6, m*4.5, [PlotDir 'Similarity_LDA_coeff_' DataSetList(nData).name], 'tif')
end

% for nData             = 1:length(DataSetList)
%     load([TempDatDir DataSetList(nData).name '.mat'])
%     numSession       = length(nDataSet3D);
%     m                = ceil(sqrt(numSession + 2));
%     figure;
%     hold on
%     nSessionData = shuffleSessionData(nDataSet, totTargets);
%     nSessionData = normalizationDim(nSessionData, 2);
%     coeffs       = coeffLDA(nSessionData +randn(size(nSessionData))*1e-3/sqrt(numTrials)* addNoise(nData), trainingTargets, testTargets);
%     imagesc(DataSetList(nData).params.timeSeries, DataSetList(nData).params.timeSeries, coeffs'*coeffs);
%     xlim([min(DataSetList(nData).params.timeSeries) max(DataSetList(nData).params.timeSeries)]);
%     ylim([min(DataSetList(nData).params.timeSeries) max(DataSetList(nData).params.timeSeries)]);
%     caxis([-1 1]);
%     axis xy;
%     gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0], 'Color','k','Linestyle','--','linewid', 0.5);
%     box off;
%     hold off;
%     xlabel('Time (s)')
%     ylabel('Time (s)')
%     setPrint(6, 4.5, [PlotDir 'Similarity_LDA_coeffAll_' DataSetList(nData).name], 'tif')
% end

%%
for nData             = 1:length(DataSetList)
    load([TempDatDir DataSetList(nData).name '.mat'])
    numSession       = length(nDataSet3D);
    m                = ceil((numSession)/4);
    numUnit          = arrayfun(@(x) length(x.nUnit), nDataSet3D, 'UniformOutput', false);
    [~, sortUnit]    = sort([numUnit{:}],'descend');    
    nDataSet3D       = nDataSet3D(sortUnit);
    figure;
    for nPlot        = 1:numSession
        subplot(m, 4, nPlot)
        hold on
        nSessionData = [permute(nDataSet3D(nPlot).unit_yes_trial,[2, 1, 3]); permute(nDataSet3D(nPlot).unit_no_trial,[2, 1, 3])];
        nSessionData = normalizationDim(nSessionData, 2);
        coeffs       = coeffPCA(nSessionData +randn(size(nSessionData))*1e-3/sqrt(numTrials)* addNoise(nData));
        imagesc(DataSetList(nData).params.timeSeries, DataSetList(nData).params.timeSeries, abs(coeffs'*coeffs));
        xlim([min(DataSetList(nData).params.timeSeries) max(DataSetList(nData).params.timeSeries)]);
        ylim([min(DataSetList(nData).params.timeSeries) max(DataSetList(nData).params.timeSeries)]);
        caxis([0 1]);
        axis xy;
        gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0], 'Color','k','Linestyle','--','linewid', 0.5);
        box off;
        hold off;
        xlabel('Time (s)')
        ylabel('Time (s)')
        title([num2str(length(nDataSet3D(nPlot).nUnit)) ' Neurons'])
    end

%     subplot(m, 4, nPlot+1)
%     caxis([0 1])
%     colorbar
%     axis off
    setPrint(4*6, m*4.5, [PlotDir 'Similarity_PCA_coeff_' DataSetList(nData).name], 'tif')
end

% for nData             = 1:length(DataSetList)
%     load([TempDatDir DataSetList(nData).name '.mat'])
%     numSession       = length(nDataSet3D);
%     m                = ceil(sqrt(numSession + 2));
%     figure;
%     hold on
%     nSessionData = shuffleSessionData(nDataSet, totTargets);
%     nSessionData = normalizationDim(nSessionData, 2);
%     coeffs       = coeffPCA(nSessionData +randn(size(nSessionData))*1e-3/sqrt(numTrials)* addNoise(nData));
%     imagesc(DataSetList(nData).params.timeSeries, DataSetList(nData).params.timeSeries, coeffs'*coeffs);
%     xlim([min(DataSetList(nData).params.timeSeries) max(DataSetList(nData).params.timeSeries)]);
%     ylim([min(DataSetList(nData).params.timeSeries) max(DataSetList(nData).params.timeSeries)]);
%     caxis([0 1]);
%     axis xy;
%     gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0], 'Color','k','Linestyle','--','linewid', 0.5);
%     box off;
%     hold off;
%     xlabel('Time (s)')
%     ylabel('Time (s)')
%     setPrint(6, 4.5, [PlotDir 'Similarity_PCA_coeffAll_' DataSetList(nData).name], 'tif')
% end

    
%%    

for nData             = 1:length(DataSetList)
    load([TempDatDir DataSetList(nData).name '.mat'])
    numSession       = length(nDataSet3D);
    m                = ceil((numSession)/4);
    numUnit          = arrayfun(@(x) length(x.nUnit), nDataSet3D, 'UniformOutput', false);
    [~, sortUnit]    = sort([numUnit{:}],'descend');    
    nDataSet3D       = nDataSet3D(sortUnit);
    figure;
    for nPlot        = 1:numSession
        subplot(m, 4, nPlot)
        hold on
        nSessionData = [permute(nDataSet3D(nPlot).unit_yes_trial,[2, 1, 3]); permute(nDataSet3D(nPlot).unit_no_trial,[2, 1, 3])];
        nSessionData = normalizationDim(nSessionData, 2);
        coeffPCAs    = coeffPCA(nSessionData +randn(size(nSessionData))*1e-3/sqrt(numTrials)* addNoise(nData));
                
        nTargets     = [true(size(nDataSet3D(nPlot).unit_yes_trial, 2),1); false(size(nDataSet3D(nPlot).unit_no_trial, 2),1)];
        nSessionData = [nSessionData; nSessionData]; %#ok<AGROW>
        coeffLDAs    = coeffLDA(nSessionData +randn(size(nSessionData))*1e-3/sqrt(numTrials)* addNoise(nData), nTargets, nTargets);

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
        title([num2str(length(nDataSet3D(nPlot).nUnit)) ' Neurons'])
    end

%     subplot(m, 4, nPlot+1)
%     caxis([0 1])
%     colorbar
%     axis off
    setPrint(4*6, m*4.5, [PlotDir 'Similarity_PCA_LDA_coeff_' DataSetList(nData).name], 'tif')
end

% for nData             = 1:length(DataSetList)
%     load([TempDatDir DataSetList(nData).name '.mat'])
%     numSession       = length(nDataSet3D);
%     m                = ceil(sqrt(numSession + 2));
%     figure;
%     hold on
%     nSessionData = shuffleSessionData(nDataSet, totTargets);
%     nSessionData = normalizationDim(nSessionData, 2);
%     coeffPCAs    = coeffPCA(nSessionData +randn(size(nSessionData))*1e-3/sqrt(numTrials)* addNoise(nData));
%     coeffLDAs    = coeffLDA(nSessionData +randn(size(nSessionData))*1e-3/sqrt(numTrials)* addNoise(nData), trainingTargets, testTargets);
%     imagesc(DataSetList(nData).params.timeSeries, DataSetList(nData).params.timeSeries, abs(coeffPCAs'*coeffLDAs));
%     xlim([min(DataSetList(nData).params.timeSeries) max(DataSetList(nData).params.timeSeries)]);
%     ylim([min(DataSetList(nData).params.timeSeries) max(DataSetList(nData).params.timeSeries)]);
%     caxis([0 1]);
%     axis xy;
%     gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0], 'Color','k','Linestyle','--','linewid', 0.5);
%     box off;
%     hold off;
%     xlabel('Time (s)')
%     ylabel('Time (s)')
%     setPrint(6, 4.5, [PlotDir 'Similarity_PCA_LDA_coeffAll_' DataSetList(nData).name], 'tif')
% end