%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 10. Similarity of PCA and LDA coefficient vectors as function of time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison_V1_10

addpath('../Func');
setDir;
load ([TempDatDir 'DataListShuffle.mat']);
addNoise         = [1 0 0 0 0 0];

if ~exist([PlotDir '/CollectedUnitsPCALDACorr'],'dir')
    mkdir([PlotDir '/CollectedUnitsPCALDACorr'])
end


numRandPickUnits    = 200;
numTrials           = numRandPickUnits*10;
totTargets          = [true(numTrials/2,1); false(numTrials/2,1)];
ROCThres            = 0.7;


% for nData             = [1 3 6] %1:length(DataSetList)
%     load([TempDatDir DataSetList(nData).name '.mat'])
%     selectedNeuronalIndex = DataSetList(nData).ActiveNeuronIndex';
%     selectedNeuronalIndex = selectedHighROCneurons(nDataSet, DataSetList(nData).params, ROCThres, selectedNeuronalIndex);
%     nDataSet              = nDataSet(selectedNeuronalIndex);
%     numUnits              = length(nDataSet);
%     figure;
%     hold on
%     nSessionData = shuffleSessionData(nDataSet(randperm(numUnits, numRandPickUnits)), totTargets, numTrials);
%     nSessionData = normalizationDim(nSessionData, 2);
%     coeffs       = coeffLDA(nSessionData, totTargets);
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
% %     setPrint(6, 4.5, [PlotDir 'CollectedUnitsPCALDACorr/Similarity_LDA_coeffAll_' DataSetList(nData).name], 'tif')
% end


for nData             = [1 3 6] %1:length(DataSetList)
    load([TempDatDir DataSetList(nData).name '.mat'])
    figure;
    hold on
    selectedNeuronalIndex = DataSetList(nData).ActiveNeuronIndex';
    selectedNeuronalIndex = selectedHighROCneurons(nDataSet, DataSetList(nData).params, ROCThres, selectedNeuronalIndex);
    nDataSet              = nDataSet(selectedNeuronalIndex);
    nSessionData = shuffleSessionData(nDataSet, totTargets, numTrials);
    nSessionData = normalizationDim(nSessionData, 2);
    coeffPCAs    = coeffPCA(nSessionData);
    coeffLDAs    = coeffLDA(nSessionData, totTargets);
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
%     setPrint(6, 4.5, [PlotDir 'CollectedUnitsPCALDACorr/Similarity_PCA_LDA_coeffAll_' DataSetList(nData).name], 'tif')
end