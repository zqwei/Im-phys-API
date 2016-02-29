%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Similarity of PCA and LDA coefficient vectors as function of time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../Func');
setDir;
load ([TempDatDir 'DataListS2CModel.mat']);
addNoise         = [1 0 0 0];

if ~exist([PlotDir 'S2CModel'],'dir')
    mkdir([PlotDir 'S2CModel'])
end

cmap                = cbrewer('div', 'Spectral', 128, 'cubic');

numFold             = 30;
numRandPickUnits    = 100;
numTrials           = numRandPickUnits*6;
totTargets          = [true(numTrials/2,1); false(numTrials/2,1)];
ROCThres            = 0.5;


for nData             = 1:length(DataSetList)
    load([TempDatDir DataSetList(nData).name '.mat'])
    selectedNeuronalIndex = DataSetList(nData).ActiveNeuronIndex';
    selectedNeuronalIndex = selectedHighROCneurons(nDataSet, DataSetList(nData).params, ROCThres, selectedNeuronalIndex);
    nDataSet              = nDataSet(selectedNeuronalIndex);
    numUnits              = length(nDataSet);
    numT                  = length(DataSetList(nData).params.timeSeries);
    corrMat               = nan(numFold, numT, numT);
    
    for nFold             = 1:numFold
        currRandPickUnits     = numRandPickUnits;
        nSessionData = shuffleSessionData(nDataSet(randperm(numUnits, currRandPickUnits)), totTargets, numTrials);
        nSessionData = normalizationDim(nSessionData, 2);
        coeffs       = coeffLDA(nSessionData, totTargets);
        corrMat(nFold, :, :) = coeffs'*coeffs;
    end
    
    
    figure;
    hold on
    
    imagesc(DataSetList(nData).params.timeSeries, DataSetList(nData).params.timeSeries, squeeze(mean(corrMat, 1)));
    xlim([min(DataSetList(nData).params.timeSeries) max(DataSetList(nData).params.timeSeries)]);
    ylim([min(DataSetList(nData).params.timeSeries) max(DataSetList(nData).params.timeSeries)]);
    caxis([0 1]);
    axis xy;
    gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0], 'Color','k','Linestyle','--','linewid', 0.5);
    box off;
    hold off;
    xlabel('LDA Time (s)')
    ylabel('LDA Time (s)')
    colormap(cmap)
    set(gca, 'TickDir', 'out')
    setPrint(8, 6, [PlotDir 'S2CModel/SimilarityLDALD_' DataSetList(nData).name])
end

close;