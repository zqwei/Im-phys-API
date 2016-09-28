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

cmap = [         0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840];


% depth
for nData                        = [1 3 4]
    if nData == 1
        load([TempDatDir DataSetList(nData).name '.mat']);
        neuronRemoveList = false(length(nDataSet), 1);
    else
        load([TempDatDir DataSetList(nData).name '_withOLRemoval.mat']);
    end
    depth                        = [DataSetList(nData).cellinfo(:).depth];    
    depth                        = depth(~neuronRemoveList);
    depthStart                   = 100;
    depthBin                     = 50;
    depthEnd                     = 800;    
    depth                        = floor((depth-depthStart)/depthBin)*depthBin+depthStart;
    depth(depth>depthEnd)        = depthEnd;
    depth(depth<depthStart)      = depthStart;
    
    uniqueDepth                  = depthStart:depthBin:depthEnd;
    
    perMat                       = nan(length(uniqueDepth), 3, numComps);
    numMat                       = nan(length(uniqueDepth), 1);
    
    for nDepth         = 1:length(uniqueDepth)
        if sum(depth == uniqueDepth(nDepth))>50  
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
                PCAmargVar(i,:)= sum((Wpca' * Xmargs{i}).^2, 2)' / totalVar;
            end
            perMat(nDepth, :, :)  = PCAmargVar(:, 1:numComps);
            numMat(nDepth)     = sum(depth == uniqueDepth(nDepth));
        end
    end
    
    figure;
    for nComps         = 1:numComps
        subplot(1, 4, nComps)
        barh(-uniqueDepth, squeeze(perMat(:, :, nComps)),'stacked', 'edgecolor', 'none')
        colormap(cmap(1:3, :))
        title(['PC' num2str(nComps)]);
        box off
        xlim([0 0.5])
        ylim([-850 0])
        set(gca, 'yTick', -800:400:0)
        xlabel('frac. EV per PC')
        ylabel('Depth (um)')
        colormap(cmap(1:3, :))
        set(gca, 'TickDir', 'out')
    end
    
    subplot(1, 4, 4)
    barh(-uniqueDepth, numMat,'k')
    title('# Units')
    box off
    set(gca, 'yTick',  -800:400:0)
    ylim([-850 0])
    xlabel('# cells')
    ylabel('Depth (um)')
    colormap(cmap(1:3, :))
    set(gca, 'TickDir', 'out')
    setPrint(8*4, 6, [PlotDir 'CollectedUnitsPCA/CollectedUnitsPCADepth_' DataSetList(nData).name '_withOLRemoval'])
end

% animal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison across animals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for nData                        = [1 4]
    if nData   == 1
        load([TempDatDir DataSetList(nData).name '.mat'])
        neuronRemoveList = false(length(nDataSet), 1);
    else
        load([TempDatDir DataSetList(nData).name '_withOLRemoval.mat'])
    end
    [~, ~, anmIndex] = unique(cell2mat({DataSetList(nData).cellinfo.anmName}'), 'rows');
    anmIndex         = anmIndex(~neuronRemoveList);
    uniAnmIndex      = unique(anmIndex);
    
    perMat                       = nan(length(uniAnmIndex), 3, numComps);
    numMat                       = nan(length(uniAnmIndex), 1);

    for nAnm    = 1:length(uniAnmIndex)
        if sum(anmIndex == uniAnmIndex(nAnm)) > 50
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
                PCAmargVar(i,:)= sum((Wpca' * Xmargs{i}).^2, 2)' / totalVar;
            end
            perMat(nAnm, :, :) = PCAmargVar(:, 1:numComps);
            numMat(nAnm)       = sum(anmIndex == uniAnmIndex(nAnm));
        end
    end

    isAnm              = ~isnan(numMat);
    perMat             = perMat(isAnm, :, :);
    numMat             = numMat(isAnm);

    for nComps         = 1:numComps
        subplot(1, 4, nComps)
        barh(1:sum(isAnm), squeeze(perMat(:, :, nComps)),'stacked', 'edgecolor', 'none')
        colormap(cmap(1:3, :))
        title(['PC' num2str(nComps)]);
        box off
        xlim([0 0.5])
        set(gca, 'yTickLabel', {})
        xlabel('frac. EV per PC')
        ylabel('Animal index')
        colormap(cmap(1:3, :))
        set(gca, 'TickDir', 'out')
    end

    subplot(1, 4, 4)
    barh(1:sum(isAnm), numMat,'k')
    title('# Units')
    box off
    set(gca, 'yTickLabel', {})
    xlabel('# cells')
    ylabel('Animal index')
    colormap(cmap(1:3, :))
    set(gca, 'TickDir', 'out')
    setPrint(8*4, 6, [PlotDir 'CollectedUnitsPCA/CollectedUnitsPCAAnm_' DataSetList(nData).name])
end

close all



load ([TempDatDir 'DataListShuffleConfounding.mat']);
nData          = 5;
load([TempDatDir DataSetList(nData).name '.mat']);
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
    PCAmargVar(i,:)= sum((Wpca' * Xmargs{i}).^2, 2)' / totalVar;
end
figure;
bar(1:numComps, PCAmargVar(:, 1:numComps)','stacked', 'edgecolor', 'none')
box off
xlim([0 numComps+0.5])
ylim([0 0.4])
xlabel('Component index')
ylabel('frac. EV per PC')
colormap(cmap(1:3, :))
set(gca, 'TickDir', 'out')
setPrint(8, 6, [PlotDir 'CollectedUnitsPCA/CollectedUnitsPCA_Pop_' DataSetList(nData).name])

APLoc          = [DataSetList(nData).cellinfo.AP_axis];
MLLoc          = [DataSetList(nData).cellinfo.ML_axis];
nDataSet       = nDataSet(APLoc>2400 & APLoc<2600 & MLLoc>1100 & MLLoc<1900);
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
    PCAmargVar(i,:)= sum((Wpca' * Xmargs{i}).^2, 2)' / totalVar;
end
figure;
bar(1:numComps, PCAmargVar(:, 1:numComps)','stacked', 'edgecolor', 'none')
box off
xlim([0 numComps+0.5])
ylim([0 0.4])
xlabel('Component index')
ylabel('frac. EV per PC')
colormap(cmap(1:3, :))
set(gca, 'TickDir', 'out')
setPrint(8, 6, [PlotDir 'CollectedUnitsPCA/CollectedUnitsPCA_Sub_' DataSetList(nData).name])
