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

ROCThres            = 0.50;


for nData      = 10 %[1 3 4 10]%1:length(DataSetList)
    if ~exist([TempDatDir DataSetList(nData).name '_withOLRemoval.mat'], 'file')
        load([TempDatDir DataSetList(nData).name '.mat'])
        neuronRemoveList = false(length(nDataSet), 1);
        selectedNeuronalIndex = DataSetList(nData).ActiveNeuronIndex';
    else
        load([TempDatDir DataSetList(nData).name '_withOLRemoval.mat'])
        selectedNeuronalIndex = DataSetList(nData).ActiveNeuronIndex(~neuronRemoveList)';
    end
    
    oldDataSet          = nDataSet;
    selectedNeuronalIndex = selectedHighROCneurons(oldDataSet, DataSetList(nData).params, ROCThres, selectedNeuronalIndex);
    disp(mean(selectedNeuronalIndex))
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
    ylim([0 0.62])
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

% % % figure;
% % % hold on
% % % for nColor = 1:length(margNames)
% % %     plot(0, nColor, 's', 'color', cmap(nColor,:), 'MarkerFaceColor',cmap(nColor,:),'MarkerSize', 8)
% % %     text(1, nColor, margNames{nColor})
% % % end
% % % xlim([0 10])
% % % hold off
% % % axis off
% % % setPrint(3, 2, [PlotDir 'CollectedUnitsPCA/CollectedUnitsPCA_Label'])
% % % 
% % % 
% % % ROCThres = 0.55;
% % % % different ROC
% % % for nData              = [1 3 4]
% % %     load([TempDatDir DataSetList(nData).name '.mat']);
% % %     selectedNeuronalIndex = DataSetList(nData).ActiveNeuronIndex';
% % %     selectedNeuronalIndex = selectedHighLocalROCneurons(nDataSet, DataSetList(nData).params, ROCThres, selectedNeuronalIndex);
% % %     sum(selectedNeuronalIndex)
% % %     nDataSet              = nDataSet(selectedNeuronalIndex);
% % %     evMat              = zeros(numFold, length(combinedParams), numComps);
% % %     firingRates        = generateDPCAData(nDataSet, numTrials);
% % %     firingRatesAverage = nanmean(firingRates, ndims(firingRates));
% % %     pcaX               = firingRatesAverage(:,:);
% % %     firingRatesAverage = bsxfun(@minus, firingRatesAverage, mean(pcaX,2));
% % %     pcaX               = bsxfun(@minus, pcaX, mean(pcaX,2));
% % %     Xmargs             = dpca_marginalize(firingRatesAverage, 'combinedParams', combinedParams, 'ifFlat', 'yes');
% % %     totalVar           = sum(sum(pcaX.^2));
% % %     [~, ~, Wpca] = svd(pcaX');
% % %     PCAmargVar         = zeros(length(combinedParams), length(nDataSet));
% % %     for i=1:length(Xmargs)
% % %         PCAmargVar(i,:) = sum((Wpca' * Xmargs{i}).^2, 2)' / totalVar;
% % %     end
% % %     
% % %     figure;
% % %     bar(1:numComps, PCAmargVar(:, 1:numComps)','stacked', 'edgecolor', 'none')
% % %     box off
% % %     xlim([0 numComps+0.5])
% % %     ylim([0 0.6])
% % %     xlabel('Component index')
% % %     ylabel('frac. EV per PC')
% % %     colormap(cmap(1:3, :))
% % %     set(gca, 'xTick', 0:5:10)
% % %     set(gca, 'TickDir', 'out')
% % %     setPrint(8, 6, [PlotDir 'CollectedUnitsPCA/CollectedUnitsPCALOCROC_' DataSetList(nData).name])
% % %     
% % % end
% % % 
% % % 
% % % 
% % % 
% % % 
% % % numComps       = 3;
% % % 
% % % 
% % % 
% % % for nData              = [1 3 4]
% % %     
% % %     load([TempDatDir DataSetList(nData).name '.mat']);
% % %     evMat              = zeros(numFold, length(combinedParams), numComps);
% % %     
% % %     figure;
% % %     
% % % %     %%% whole population
% % % %     firingRates        = generateDPCAData(nDataSet, numTrials);
% % % %     firingRatesAverage = nanmean(firingRates, ndims(firingRates));
% % % % 
% % % %     pcaX               = firingRatesAverage(:,:);
% % % %     firingRatesAverage = bsxfun(@minus, firingRatesAverage, mean(pcaX,2));
% % % %     pcaX               = bsxfun(@minus, pcaX, mean(pcaX,2));
% % % % 
% % % %     % marginalizing
% % % %     Xmargs             = dpca_marginalize(firingRatesAverage, 'combinedParams', combinedParams, 'ifFlat', 'yes');
% % % %     % total variance
% % % %     totalVar           = sum(sum(pcaX.^2));
% % % % 
% % % %     % PCA explained variance
% % % %     [~, ~, Wpca] = svd(pcaX');
% % % % 
% % % %     PCAmargVar         = zeros(length(combinedParams), length(nDataSet));
% % % % 
% % % %     for i=1:length(Xmargs)
% % % %         PCAmargVar(i,:) = sum((Wpca' * Xmargs{i}).^2, 2)' / totalVar * 100;
% % % %     end
% % % %     PCAmargVar          = bsxfun(@rdivide, PCAmargVar, sum(PCAmargVar));
% % % % 
% % % %     subplot(1, 4, 1);
% % % % %     plot(1:numComps, PCAmargVar(:, 1:numComps)', '-o', 'linewid', 2)
% % % %     bar(1:numComps, PCAmargVar(:, 1:numComps)','stacked', 'edgecolor', 'none')
% % % %     box off
% % % %     xlim([0.5 3.5])
% % % %     ylim([0 1])
% % % %     legend(margNames)
% % % % %     legend('location','best')
% % % %     legend boxoff
% % % %     xlabel('Component index')
% % % %     ylabel('frac. EV per PC')
% % % %     title('Whole population')
% % % %     set(gca, 'TickDir', 'out')
% % %     
% % %     
% % %     numSubUnit         = 100;
% % %     for nFold          = 1:numFold        
% % %         randIndex          = randperm(length(nDataSet), numSubUnit);    
% % %         firingRates        = generateDPCAData(nDataSet(randIndex), numTrials);
% % %         firingRatesAverage = nanmean(firingRates, ndims(firingRates));
% % % 
% % %         pcaX               = firingRatesAverage(:,:);
% % %         firingRatesAverage = bsxfun(@minus, firingRatesAverage, mean(pcaX,2));
% % %         pcaX               = bsxfun(@minus, pcaX, mean(pcaX,2));
% % % 
% % %         % marginalizing
% % %         Xmargs             = dpca_marginalize(firingRatesAverage, 'combinedParams', combinedParams, 'ifFlat', 'yes');
% % %         % total variance
% % %         totalVar           = sum(sum(pcaX.^2));
% % % 
% % %         % PCA explained variance
% % %         [~, ~, Wpca] = svd(pcaX');
% % % 
% % %         PCAmargVar         = zeros(length(combinedParams), numSubUnit);
% % % 
% % %         for i=1:length(Xmargs)
% % %             PCAmargVar(i,:) = sum((Wpca' * Xmargs{i}).^2, 2)' / totalVar * 100;
% % %         end
% % %         PCAmargVar          = bsxfun(@rdivide, PCAmargVar, sum(PCAmargVar));
% % % 
% % %         evMat(nFold, :, :)  = PCAmargVar(:, 1:numComps);
% % %     end
% % % 
% % %     subplot(1, 4, 2);
% % %     errorbar((1:numComps)'*ones(length(combinedParams), 1)', squeeze(mean(evMat))', squeeze(std(evMat))'/sqrt(numFold), '-o', 'linewid', 2)
% % %     box off
% % %     xlim([0.5 3.5])
% % %     ylim([0 1])
% % % %     legend(margNames)
% % %     xlabel('Component index')
% % %     ylabel('frac. EV per PC')
% % %     title(['Subpopulation n=' num2str(numSubUnit)])
% % %     set(gca, 'TickDir', 'out')
% % %     
% % %     numSubUnit         = 200;
% % %     for nFold          = 1:numFold        
% % %         randIndex          = randperm(length(nDataSet), numSubUnit);    
% % %         firingRates        = generateDPCAData(nDataSet(randIndex), numTrials);
% % %         firingRatesAverage = nanmean(firingRates, ndims(firingRates));
% % % 
% % %         pcaX               = firingRatesAverage(:,:);
% % %         firingRatesAverage = bsxfun(@minus, firingRatesAverage, mean(pcaX,2));
% % %         pcaX               = bsxfun(@minus, pcaX, mean(pcaX,2));
% % % 
% % %         % marginalizing
% % %         Xmargs             = dpca_marginalize(firingRatesAverage, 'combinedParams', combinedParams, 'ifFlat', 'yes');
% % %         % total variance
% % %         totalVar           = sum(sum(pcaX.^2));
% % % 
% % %         % PCA explained variance
% % %         [~, ~, Wpca] = svd(pcaX');
% % % 
% % %         PCAmargVar         = zeros(length(combinedParams), numSubUnit);
% % % 
% % %         for i=1:length(Xmargs)
% % %             PCAmargVar(i,:) = sum((Wpca' * Xmargs{i}).^2, 2)' / totalVar * 100;
% % %         end
% % %         PCAmargVar          = bsxfun(@rdivide, PCAmargVar, sum(PCAmargVar));
% % % 
% % %         evMat(nFold, :, :)  = PCAmargVar(:, 1:numComps);
% % %     end
% % % 
% % %     subplot(1, 4, 3);
% % %     errorbar((1:numComps)'*ones(length(combinedParams), 1)', squeeze(mean(evMat))', squeeze(std(evMat))'/sqrt(numFold), '-o', 'linewid', 2)
% % %     box off
% % %     xlim([0.5 3.5])
% % %     ylim([0 1])
% % % %     legend(margNames)
% % %     xlabel('Component index')
% % %     ylabel('frac. EV per PC')
% % %     title(['Subpopulation n=' num2str(numSubUnit)])
% % %     set(gca, 'TickDir', 'out')
% % %     
% % %     numSubUnit         = 500;
% % %     for nFold          = 1:numFold        
% % %         randIndex          = randperm(length(nDataSet), numSubUnit);    
% % %         firingRates        = generateDPCAData(nDataSet(randIndex), numTrials);
% % %         firingRatesAverage = nanmean(firingRates, ndims(firingRates));
% % % 
% % %         pcaX               = firingRatesAverage(:,:);
% % %         firingRatesAverage = bsxfun(@minus, firingRatesAverage, mean(pcaX,2));
% % %         pcaX               = bsxfun(@minus, pcaX, mean(pcaX,2));
% % % 
% % %         % marginalizing
% % %         Xmargs             = dpca_marginalize(firingRatesAverage, 'combinedParams', combinedParams, 'ifFlat', 'yes');
% % %         % total variance
% % %         totalVar           = sum(sum(pcaX.^2));
% % % 
% % %         % PCA explained variance
% % %         [~, ~, Wpca] = svd(pcaX');
% % % 
% % %         PCAmargVar         = zeros(length(combinedParams), numSubUnit);
% % % 
% % %         for i=1:length(Xmargs)
% % %             PCAmargVar(i,:) = sum((Wpca' * Xmargs{i}).^2, 2)' / totalVar * 100;
% % %         end
% % %         PCAmargVar          = bsxfun(@rdivide, PCAmargVar, sum(PCAmargVar));
% % % 
% % %         evMat(nFold, :, :)  = PCAmargVar(:, 1:numComps);
% % %     end
% % % 
% % %     subplot(1, 4, 4);
% % %     errorbar((1:numComps)'*ones(length(combinedParams), 1)', squeeze(mean(evMat))', squeeze(std(evMat))'/sqrt(numFold), '-o', 'linewid', 2)
% % %     box off
% % %     xlim([0.5 3.5])
% % %     ylim([0 1])
% % % %     legend(margNames)
% % %     xlabel('Component index')
% % %     ylabel('frac. EV per PC')
% % %     title(['Subpopulation n=' num2str(numSubUnit)])
% % %     set(gca, 'TickDir', 'out')
% % %     
% % %     setPrint(8*4, 6, [PlotDir 'CollectedUnitsPCA/CollectedSubUnitsPCA_' DataSetList(nData).name])
% % %     
% % % end


close all