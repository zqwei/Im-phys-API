%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Similarity of PCA and LDA coefficient vectors as function of time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../Func');
setDir;
load ([TempDatDir 'DataListShuffle.mat']);
addNoise         = [1 0 0 0];

if ~exist([PlotDir 'CollectedUnitsPCALDACorr'],'dir')
    mkdir([PlotDir 'CollectedUnitsPCALDACorr'])
end


numRandPickUnits    = 100;
numTrials           = numRandPickUnits*10;
totTargets          = [true(numTrials/2,1); false(numTrials/2,1)];
ROCThres            = 0.5;


for nData             = 1:length(DataSetList)-1
    load([TempDatDir DataSetList(nData).name '.mat'])
    selectedNeuronalIndex = DataSetList(nData).ActiveNeuronIndex';
    selectedNeuronalIndex = selectedHighROCneurons(nDataSet, DataSetList(nData).params, ROCThres, selectedNeuronalIndex);
    nDataSet              = nDataSet(selectedNeuronalIndex);
    numUnits              = length(nDataSet);
    figure;
    hold on
    currRandPickUnits     = numRandPickUnits;
    if currRandPickUnits>numUnits; currRandPickUnits = 100; end
    nSessionData = shuffleSessionData(nDataSet(randperm(numUnits, currRandPickUnits)), totTargets, numTrials);
    nSessionData = normalizationDim(nSessionData, 2);
    coeffs       = coeffLDA(nSessionData, totTargets);
    imagesc(DataSetList(nData).params.timeSeries, DataSetList(nData).params.timeSeries, coeffs'*coeffs);
    xlim([min(DataSetList(nData).params.timeSeries) max(DataSetList(nData).params.timeSeries)]);
    ylim([min(DataSetList(nData).params.timeSeries) max(DataSetList(nData).params.timeSeries)]);
%     caxis([-1 1]);
    colorbar
    axis xy;
    gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0], 'Color','k','Linestyle','--','linewid', 0.5);
    box off;
    hold off;
    xlabel('Time (s)')
    ylabel('Time (s)')
    setPrint(8, 6, [PlotDir 'CollectedUnitsPCALDACorr/SimilarityLDALDA_' DataSetList(nData).name], 'pdf')
end
