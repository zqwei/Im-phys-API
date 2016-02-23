%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dPCA and PCA across trials
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../Func');
setDir;
load ([TempDatDir 'DataListShuffle.mat']);

combinedParams = {{1}, {2}, {[1 2]}};
margNames      = {'Stim', 'Time', 'Inter'};
numTrials      = 100;


numComps             = 15;
if ~exist([PlotDir 'CollectedUnitsdPCA'],'dir')
    mkdir([PlotDir 'CollectedUnitsdPCA'])
end

cmap = cbrewer('qual', 'Set1', 3, 'cubic');

for nData              = 1:length(DataSetList)-1
    load([TempDatDir DataSetList(nData).name '.mat']);
    time               = DataSetList(nData).params.timeSeries;
    timeEvents         = [DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0];
    firingRates        = generateDPCAData(nDataSet, numTrials);
    firingRatesAverage = nanmean(firingRates, ndims(firingRates));
    trialNum           = ones(size(firingRatesAverage, 1), size(firingRatesAverage, 2))*numTrials;
    optimalLambda = dpca_optimizeLambda(firingRatesAverage, firingRates, ...
                trialNum, ...
                'numComps', numComps, ....
                'combinedParams', combinedParams, ...
                'numRep', 10, ...  % increase this number to ~10 for better accuracy
                'filename', [TempDatDir 'optimalLambdas_' DataSetList(nData).name '.mat'],...
                'display','no');
    [W,V,whichMarg] = dpca(firingRatesAverage, numComps, ...
                'combinedParams', combinedParams, ...
                'lambda', optimalLambda);
    explVar = dpcaEV(firingRatesAverage, W, 'combinedParams', combinedParams);
    
    save([TempDatDir 'dPCA_' DataSetList(nData).name '.mat'], 'optimalLambda', 'W', 'V', 'explVar')
    
    load([TempDatDir 'dPCA_' DataSetList(nData).name '.mat'], 'explVar')
    
%     figure;
%     subplot(1, 3, 1)
%     plot(1:numComps, [explVar.PCAcomponentVar(1:numComps)', explVar.dPCAcomponentVar'], '-o');
%     legend({'PCA', 'dPCA'})
%     box off
%     xlim([0 numComps+1])
%     set(gca,'xTick', 0:5:numComps)
%     xlabel('Component index')
%     ylabel('% EV')
%     
%     subplot(1, 3, 2)
%     bar(1:numComps, explVar.PCAmargVar(:, 1:numComps)', 'stack', 'edgecolor', 'none');
%     box off
%     legend(margNames)
%     xlim([0 numComps+1])
%     xlabel('Component index')
%     set(gca,'xTick', 0:5:numComps)
%     ylabel('% EV per PC')
%     
%     subplot(1, 3, 3)
%     bar(1:numComps, explVar.dPCAmargVar', 'stack', 'edgecolor', 'none');
%     box off
%     legend(margNames)
%     xlim([0 numComps+1])
%     xlabel('Component index')
%     set(gca,'xTick', 0:5:numComps)
%     ylabel('% EV per dPC')
%     
%     setPrint(8*3, 6, [PlotDir 'CollectedUnitsdPCA/CollectedUnitsdPCA_' DataSetList(nData).name])
    
    figure
    bar(1:numComps, explVar.PCAmargVar(:, 1:numComps)', 'stack', 'edgecolor', 'none');
    box off
    legend(margNames)
    xlim([0 numComps+1])
    xlabel('Component index')
    set(gca,'xTick', 0:5:numComps)
    ylabel('% EV per PC')
    colormap(cmap)
    setPrint(8, 6, [PlotDir 'CollectedUnitsdPCA/CollectedUnitsPCA_' DataSetList(nData).name])
    
end

close all