%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dPCA and PCA across trials
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../Func');
setDir;
load ([TempDatDir 'DataListShuffle.mat']);

combinedParams = {{1}, {2}, {[1 2]}};
margNames      = {'Stim', 'Time', 'Inter'};
numTrials      = 100;
numComps       = 3;

if ~exist([PlotDir 'CollectedUnitsPCA'],'dir')
    mkdir([PlotDir 'CollectedUnitsPCA'])
end

for nData              = [1 3 4]
    if nData   == 1
        load([TempDatDir DataSetList(nData).name '.mat'])
        neuronRemoveList = false(length(nDataSet), 1);
    else
        load([TempDatDir DataSetList(nData).name '_withOLRemoval.mat'])
    end
    firingRates        = generateDPCAData(nDataSet, numTrials);
    firingRatesAverage = nanmean(firingRates, ndims(firingRates));
    pcaFiringRatesAverage = zeros(numComps, 2, 77);
    firingRatesAverage = [squeeze(firingRatesAverage(:, 1, :)), squeeze(firingRatesAverage(:, 2, :));];
    [~,score,~]        = pca(firingRatesAverage', 'NumComponents', numComps);
    pcaFiringRatesAverage(:, 1, :) = score(1:77, :)';
    pcaFiringRatesAverage(:, 2, :) = score(78:end, :)';
    
    figure;
    hold on
    plot(squeeze(pcaFiringRatesAverage(1, 1, :)), squeeze(pcaFiringRatesAverage(2, 1, :)), '-b', 'linewid', 2);
    plot(squeeze(pcaFiringRatesAverage(1, 2, :)), squeeze(pcaFiringRatesAverage(2, 2, :)),  '-r', 'linewid', 2);
    hold off
    box off
    xlabel('PC1 score')
    ylabel('PC2 score')  
    set(gca, 'TickDir', 'out')
    
    setPrint(8, 6, [PlotDir 'CollectedUnitsPCA/CollectedUnitsPCATrace2D_' DataSetList(nData).name])
    
end
