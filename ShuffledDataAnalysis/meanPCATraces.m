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

if ~exist([PlotDir 'CollectedUnitsdPCA'],'dir')
    mkdir([PlotDir 'CollectedUnitsdPCA'])
end

for nData              = 1:length(DataSetList)-1
    load([TempDatDir DataSetList(nData).name '.mat']);
    firingRates        = generateDPCAData(nDataSet, numTrials);
    firingRatesAverage = nanmean(firingRates, ndims(firingRates));
    pcaFiringRatesAverage = zeros(numComps, 2, 77);
    firingRatesAverage = [squeeze(firingRatesAverage(:, 1, :)), squeeze(firingRatesAverage(:, 2, :));];
    [~,score,~]        = pca(firingRatesAverage', 'NumComponents', numComps);
    pcaFiringRatesAverage(:, 1, :) = score(1:77, :)';
    pcaFiringRatesAverage(:, 2, :) = score(78:end, :)';
    
    figure;
    for nPlot           = 1:numComps
        subplot(1, numComps, nPlot)
        hold on
        plot(DataSetList(nData).params.timeSeries, squeeze(pcaFiringRatesAverage(nPlot, :, :)), '-', 'linewid', 2);
        gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 0.5)
        hold off
        box off
        xlim([min(DataSetList(nData).params.timeSeries) max(DataSetList(nData).params.timeSeries)]);
        set(gca,'xTick', 0:5:numComps)
        xlabel('Time (s)')
        ylabel('PC score')  
    end
    
    setPrint(8*3, 6, [PlotDir 'CollectedUnitsdPCA/CollectedUnitsdPCATrace_' DataSetList(nData).name], 'pdf')
    
end

close all