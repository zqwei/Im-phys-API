%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Similarity of PCA and LDA coefficient vectors as function of time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../Func');
setDir;
load ([TempDatDir 'DataListModeled.mat']);

if ~exist([PlotDir 'ModeledCollectedUnitsPCALDACorr'],'dir')
    mkdir([PlotDir 'ModeledCollectedUnitsPCALDACorr'])
end


numRandPickUnits    = 400;
numTrials           = numRandPickUnits*3;
totTargets          = [true(numTrials/2,1); false(numTrials/2,1)];
ROCThres            = 0.5;


for nData             =1:length(DataSetList)
    load([TempDatDir DataSetList(nData).name '.mat'])
    figure;
    hold on
    selectedNeuronalIndex = DataSetList(nData).ActiveNeuronIndex';
    selectedNeuronalIndex = selectedHighROCneurons(nDataSet, DataSetList(nData).params, ROCThres, selectedNeuronalIndex);
    nDataSet              = nDataSet(selectedNeuronalIndex);
    numUnits              = length(nDataSet);
    currRandPickUnits     = numRandPickUnits;
    if currRandPickUnits>numUnits; currRandPickUnits = 100; end
    nSessionData = shuffleSessionData(nDataSet(randperm(numUnits, currRandPickUnits)), totTargets, numTrials);
%     nSessionData = normalizationDim(nSessionData, 2);
    coeffPCAs    = coeffPCA(nSessionData);
    coeffLDAs    = coeffLDA(nSessionData, totTargets);
    imagesc(DataSetList(nData).params.timeSeries, DataSetList(nData).params.timeSeries, abs(coeffPCAs'*coeffLDAs));
    xlim([min(DataSetList(nData).params.timeSeries) max(DataSetList(nData).params.timeSeries)]);
    ylim([min(DataSetList(nData).params.timeSeries) max(DataSetList(nData).params.timeSeries)]);
%     caxis([0 0.6]);
    colorbar
    axis xy;
    gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0], 'Color','k','Linestyle','--','linewid', 0.5);
    box off;
    hold off;
    xlabel('Time (s)')
    ylabel('Time (s)')
    setPrint(8, 6, [PlotDir 'CollectedUnitsPCALDACorr/SimilarityLDAPCANONormalization_' DataSetList(nData).name], 'pdf')
end
