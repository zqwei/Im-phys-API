%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dPCA and PCA across trials
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

ROCThres            = 0.70;


for nData      = [1 10]
    if ~exist([TempDatDir DataSetList(nData).name '_withOLRemoval.mat'], 'file')
        load([TempDatDir DataSetList(nData).name '.mat'])
        neuronRemoveList = false(length(nDataSet), 1);
        selectedNeuronalIndex = DataSetList(nData).ActiveNeuronIndex';
    else
        load([TempDatDir DataSetList(nData).name '_withOLRemoval.mat'])
        selectedNeuronalIndex = DataSetList(nData).ActiveNeuronIndex(~neuronRemoveList)';
    end
    
    depth_list          = [nDataSet.depth_in_um]';
    
    oldDataSet          = nDataSet;
    selectedNeuronalIndex = selectedHighROCneurons(oldDataSet, DataSetList(nData).params, ROCThres, selectedNeuronalIndex);
    selectedNeuronalIndex = selectedNeuronalIndex & depth_list < 471;
    nDataSet              = oldDataSet(selectedNeuronalIndex);
    decodability          = zeros(numFold, size(nDataSet(1).unit_yes_trial,2)); 
    
    evMat              = zeros(numFold, length(combinedParams), numComps);
    firingRates        = generateDPCAData(nDataSet, numTrials);
    firingRatesAverage = nanmean(firingRates, ndims(firingRates));
    pcaX               = firingRatesAverage(:,:);
    firingRatesAverage = bsxfun(@minus, firingRatesAverage, mean(pcaX,2));
    pcaX               = bsxfun(@minus, pcaX, mean(pcaX,2));
    Xmargs             = dpca_marginalize(firingRatesAverage, 'combinedParams', combinedParams, 'ifFlat', 'yes');
    totalVar           = sum(sum(pcaX.^2));
    [~, S, Wpca] = svd(pcaX');

    PCAmargVar         = zeros(length(combinedParams), length(nDataSet));
    for i=1:length(Xmargs)
        PCAmargVar(i,:) = sum((Wpca' * Xmargs{i}).^2, 2)' / totalVar;
    end
    sum(sum(PCAmargVar, 1)>0.01)

    figure;
    bar(1:numComps, PCAmargVar(:, 1:numComps)','stacked', 'edgecolor', 'none')
    box off
    xlim([0 numComps+0.5])
    ylim([0 0.5])
    xlabel('Component index')
    ylabel('frac. EV per PC')
    colormap(cmap(1:3, :))
    set(gca, 'xTick', 0:5:10)
    set(gca, 'TickDir', 'out')
    setPrint(8, 6, [PlotDir 'CollectedUnitsPCA/CollectedUnitsPCA_' DataSetList(nData).name])  
    legend({'Trial type', 'Time', 'Other'})
    legend('boxoff')
    setPrint(8, 6, [PlotDir 'CollectedUnitsPCA/' DataSetList(nData).name '_pca'], 'svg')   
end



for nData                        = 10
    if nData == 1
        load([TempDatDir DataSetList(nData).name '.mat']);
        neuronRemoveList = false(length(nDataSet), 1);
    else
        load([TempDatDir DataSetList(nData).name '_withOLRemoval.mat']);
    end
    
    depth_list          = [nDataSet.depth_in_um]';
    
    oldDataSet          = nDataSet;
    selectedNeuronalIndex = selectedHighROCneurons(oldDataSet, DataSetList(nData).params, ROCThres, selectedNeuronalIndex);
    selectedNeuronalIndex = selectedNeuronalIndex & depth_list < 471;
    nDataSet              = oldDataSet(selectedNeuronalIndex);
    
    depth                        = [nDataSet.depth_in_um];    
    depthStart                   = 100;
    depthBin                     = 100;
    depthEnd                     = 800;    
    depth                        = floor((depth-depthStart)/depthBin)*depthBin+depthStart;
    depth(depth>depthEnd)        = depthEnd;
    depth(depth<depthStart)      = depthStart;
    
    uniqueDepth                  = depthStart:depthBin:depthEnd;
    
    perMat                       = nan(length(uniqueDepth), 3, numComps);
    numMat                       = nan(length(uniqueDepth), 1);
    
    for nDepth         = 1:length(uniqueDepth)
        if sum(depth == uniqueDepth(nDepth))>10  
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
        xlim([0 0.9])
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