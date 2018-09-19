addpath('../Func');
setDir;
load ([TempDatDir 'DataListShuffle.mat']);

combinedParams = {{1}, {2}, {[1 2]}};
margNames      = {'Stim', 'Time', 'Inter'};
numTrials      = 100;
numFold        = 100;
if ~exist([PlotDir 'CollectedUnitsPCA'],'dir')
    mkdir([PlotDir 'CollectedUnitsPCA'])
end


numComps       = 10;
cmap = [         0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840];


for nData              = [1 3 4]
    if nData   == 1
        load([TempDatDir DataSetList(nData).name '.mat'])
        neuronRemoveList = false(length(nDataSet), 1);
    else
        load([TempDatDir DataSetList(nData).name '_withOLRemoval.mat'])
    end

    evMat              = zeros(numFold, length(combinedParams), numComps);
    firingRates        = generateDPCAData(nDataSet, numTrials);
    firingRatesAverage = nanmean(firingRates, ndims(firingRates));
    pcaX               = firingRatesAverage(:,:);
    firingRatesAverage = bsxfun(@minus, firingRatesAverage, mean(pcaX,2));
    pcaX               = bsxfun(@minus, pcaX, mean(pcaX,2));
    Xmargs             = dpca_marginalize(firingRatesAverage, 'combinedParams', combinedParams, 'ifFlat', 'yes');
    totalVar           = sum(sum(pcaX.^2));
    [~, ~, Wpca] = svd(pcaX');
    X = pcaX';
    W = Wpca';
    T = X * W(:, 1:6);
    pcX = T * W(:, 1:6)';
    
    [bestEV, bestFit] = min(mean((pcX - X).^2, 1)./var(X, [], 1)); %
    % [bestEV, bestFit] = max(diag(corr(pcX, X, 'type', 'Spearman')));
    disp(bestEV)
    firingRatesAverage = nanmean(firingRates, ndims(firingRates));
    pcaX               = firingRatesAverage(:,:);
    aveFit = mean(pcaX, 2);
    
    figure;
    subplot(1, 2, 1)
    hold on
    plot(DataSetList(nData).params.timeSeries, pcaX(bestFit, 1:2:end), '-b', 'linewid', 2);
    plot(DataSetList(nData).params.timeSeries, pcaX(bestFit, 2:2:end), '-r', 'linewid', 2);
%     ylim([0 150])
    gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 0.5)
    hold off
    box off
    xlim([DataSetList(nData).params.timeSeries(2) DataSetList(nData).params.timeSeries(end-1)]);
    xlabel('Time (s)')
    ylabel(['Spikes /s'])  
    set(gca, 'TickDir', 'out')
    
    subplot(1, 2, 2)
    hold on
    plot(DataSetList(nData).params.timeSeries, pcX(1:2:end, bestFit) + aveFit(bestFit), '-b', 'linewid', 2);
    plot(DataSetList(nData).params.timeSeries, pcX(2:2:end, bestFit) + aveFit(bestFit), '-r', 'linewid', 2);
%     ylim([0 150])
    gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 0.5)
    hold off
    box off
    xlim([DataSetList(nData).params.timeSeries(2) DataSetList(nData).params.timeSeries(end-1)]);
    xlabel('Time (s)')
    ylabel(['Spikes /s'])  
    set(gca, 'TickDir', 'out')

end