%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dPCA and PCA across trials
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../Func');
setDir;
load ([TempDatDir 'DataListEphys.mat']);

combinedParams = {{1}, {2}, {[1 2]}};
margNames      = {'Stim', 'Time', 'Inter'};
numTrials      = 100;
numFold        = 100;


numComps       = 10;
cmap = [         0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840];


for nData              = 1:length(DataSetList)
    load([TempDatDir DataSetList(nData).name '.mat']);
    evMat              = zeros(numFold, length(combinedParams), numComps);
    firingRates        = generateDPCAData(nDataSet, numTrials);
    firingRatesAverage = nanmean(firingRates, ndims(firingRates));
    pcaX               = firingRatesAverage(:,:);
    firingRatesAverage = bsxfun(@minus, firingRatesAverage, mean(pcaX,2));
    pcaX               = bsxfun(@minus, pcaX, mean(pcaX,2));
    Xmargs             = dpca_marginalize(firingRatesAverage, 'combinedParams', combinedParams, 'ifFlat', 'yes');
    totalVar           = sum(sum(pcaX.^2));
    [~, ~, Wpca] = svd(pcaX');
    PCAmargVar         = zeros(length(combinedParams), length(nDataSet));
    for i=1:length(Xmargs)
        PCAmargVar(i,:) = sum((Wpca' * Xmargs{i}).^2, 2)' / totalVar;
    end
    
    figure;
    bar(1:numComps, PCAmargVar(:, 1:numComps)','stacked')
    box off
    xlim([0 numComps+0.5])
    ylim([0 0.5])
    xlabel('Component index')
    ylabel('frac. EV per PC')
    colormap(cmap(1:3, :))
    set(gca, 'xTick', 0:5:10)
    set(gca, 'TickDir', 'out')
    setPrint(8, 6, [PlotDir 'EphysBinningTest/CollectedUnitsPCA_' DataSetList(nData).name])
    
end


close all