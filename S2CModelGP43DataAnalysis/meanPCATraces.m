%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dPCA and PCA across trials
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../Func');
setDir;
load ([TempDatDir 'DataListS2CGP43Model.mat']);

combinedParams = {{1}, {2}, {[1 2]}};
margNames      = {'Stim', 'Time', 'Inter'};
numTrials      = 100;
numComps       = 3;


for nData              = 2
    load([TempDatDir DataSetList(nData).name '.mat']);
    firingRates        = generateDPCAData(nDataSet, numTrials);
    firingRatesAverage = nanmean(firingRates, ndims(firingRates));
    pcaFiringRatesAverage = zeros(numComps, 2, 77);
    firingRatesAverage = [squeeze(firingRatesAverage(:, 1, :)), squeeze(firingRatesAverage(:, 2, :));];
    [coeffs,score,latent]        = pca(firingRatesAverage', 'NumComponents', numComps);
    pcaFiringRatesAverage(:, 1, :) = score(1:77, :)';
    pcaFiringRatesAverage(:, 2, :) = score(78:end, :)';
    
    figure;
    for nPlot           = 1:numComps
        subplot(1, numComps, nPlot)
        hold on
        plot(DataSetList(nData).params.timeSeries, squeeze(pcaFiringRatesAverage(nPlot, 1, :)), '-r', 'linewid', 2);
        plot(DataSetList(nData).params.timeSeries, squeeze(pcaFiringRatesAverage(nPlot, 2, :)), '-b', 'linewid', 2);
        gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 0.5)
        hold off
        box off
        xlim([min(DataSetList(nData).params.timeSeries) max(DataSetList(nData).params.timeSeries)]);
        xlabel('Time (s)')
        ylabel(['PC' num2str(nPlot) ' score'])  
    end
    
    setPrint(8*3, 6, [PlotDir 'S2CGP43Model/CollectedUnitsPCATrace_' DataSetList(nData).name])
    
end


close all