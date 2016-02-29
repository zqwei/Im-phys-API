%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Collected population PCA variance and EV of trial type over time
%
% fixed number of neuron with normalization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../Func');
setDir;
load ([TempDatDir 'DataListShuffle.mat']);

numUnits       = 100;
numTrials      = numUnits * 3;
numComps       = 2;
ROCThres       = 0.5;
trialType      = [true(numTrials, 1); false(numTrials, 1)];
numFold        = 30;

if ~exist([PlotDir 'CollectedUnitsPCA'],'dir')
    mkdir([PlotDir 'CollectedUnitsPCA'])
end



for nData              = [1 3 4]
    load([TempDatDir DataSetList(nData).name '.mat']);
    numT               = size(nDataSet(1).unit_yes_trial, 2);
    pcaVar             = nan(numFold, numComps, numT);
    evTrialType        = nan(numFold, numComps, numT);
    selectedNeuronalIndex = DataSetList(nData).ActiveNeuronIndex';
    selectedNeuronalIndex = selectedHighROCneurons(nDataSet, DataSetList(nData).params, ROCThres, selectedNeuronalIndex);
    nDataSet              = nDataSet(selectedNeuronalIndex);

    for nFold = 1:numFold
        currNumUnits = numUnits;
        randPickUnits         = randperm(length(nDataSet));
        randPickUnits         = randPickUnits(1:currNumUnits);
        firingRates           = generatePCAData(nDataSet(randPickUnits), numTrials);

        for nTime          = 1:numT
            nFiringRates       = squeeze(firingRates(nTime, :, :))';
            nFiringRates       = bsxfun(@minus, nFiringRates, mean(nFiringRates));
            nFiringRates       = bsxfun(@rdivide, nFiringRates, std(nFiringRates));
            nFiringRates(isnan(nFiringRates)) = 0;
            [~,score,latent,~,explained,~] = pca(nFiringRates, 'NumComponents', numComps);
            totVar            = sum(latent);
            pcaVar(nFold, :, nTime)  = explained(1:numComps)/sum(explained);
            evTrialType(nFold, :, nTime) = mean(score(1:numTrials, :)).^2/totVar;
        end        
    end
        
    figure;
    subplot(2, 1, 1)
    hold on
    
    shadedErrorBar(DataSetList(nData).params.timeSeries, squeeze(mean(pcaVar(:, 1, :), 1)),...
            squeeze(std(pcaVar(:, 1, :), [], 1))/sqrt(numFold),...
            {'-', 'linewid', 1.0, 'color', 'k'}, 0.5);
        
    shadedErrorBar(DataSetList(nData).params.timeSeries, squeeze(mean(pcaVar(:, 2, :), 1)),...
            squeeze(std(pcaVar(:, 2, :), [], 1))/sqrt(numFold),...
            {'-', 'linewid', 1.0, 'color', 'm'}, 0.5);  
    legend({'PC1', 'PC2'})    
    axis xy;
    xlim([min(DataSetList(nData).params.timeSeries) max(DataSetList(nData).params.timeSeries)]);
    ylim([0 0.5])
    gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 0.5)
    box off;
    hold off;
    xlabel('Time (s)');
    ylabel('frac. PC var.');
    set(gca, 'TickDir', 'out')

    subplot(2, 1, 2)
    hold on
    shadedErrorBar(DataSetList(nData).params.timeSeries, squeeze(mean(evTrialType(:, 1, :), 1)),...
            squeeze(std(evTrialType(:, 1, :), [], 1))/sqrt(numFold),...
            {'-', 'linewid', 1.0, 'color', 'k'}, 0.5);
        
    shadedErrorBar(DataSetList(nData).params.timeSeries, squeeze(mean(evTrialType(:, 2, :), 1)),...
            squeeze(std(evTrialType(:, 2, :), [], 1))/sqrt(numFold),...
            {'-', 'linewid', 1.0, 'color', 'm'}, 0.5);  
    legend({'PC1', 'PC2'})    
    axis xy;
    xlim([min(DataSetList(nData).params.timeSeries) max(DataSetList(nData).params.timeSeries)]);
    ylim([0 0.5])
    gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 0.5)
    box off;
    hold off;
    xlabel('Time (s)');
    ylabel('frac. EV stim.');
    set(gca, 'TickDir', 'out')
    
    setPrint(2, 6*2, [PlotDir 'CollectedUnitsPCALDACorr/SimilarityPCA_100_' DataSetList(nData).name])
end


numUnits       = 500;
numTrials      = numUnits * 3;
numComps       = 2;
ROCThres       = 0.5;
trialType      = [true(numTrials, 1); false(numTrials, 1)];
numFold        = 30;

if ~exist([PlotDir 'CollectedUnitsPCA'],'dir')
    mkdir([PlotDir 'CollectedUnitsPCA'])
end



for nData              = [1 3 4]
    load([TempDatDir DataSetList(nData).name '.mat']);
    numT               = size(nDataSet(1).unit_yes_trial, 2);
    pcaVar             = nan(numFold, numComps, numT);
    evTrialType        = nan(numFold, numComps, numT);
    selectedNeuronalIndex = DataSetList(nData).ActiveNeuronIndex';
    selectedNeuronalIndex = selectedHighROCneurons(nDataSet, DataSetList(nData).params, ROCThres, selectedNeuronalIndex);
    nDataSet              = nDataSet(selectedNeuronalIndex);

    for nFold = 1:numFold
        currNumUnits = numUnits;
        randPickUnits         = randperm(length(nDataSet));
        randPickUnits         = randPickUnits(1:currNumUnits);
        firingRates           = generatePCAData(nDataSet(randPickUnits), numTrials);

        for nTime          = 1:numT
            nFiringRates       = squeeze(firingRates(nTime, :, :))';
            nFiringRates       = bsxfun(@minus, nFiringRates, mean(nFiringRates));
            nFiringRates       = bsxfun(@rdivide, nFiringRates, std(nFiringRates));
            nFiringRates(isnan(nFiringRates)) = 0;
            [~,score,latent,~,explained,~] = pca(nFiringRates, 'NumComponents', numComps);
            totVar            = sum(latent);
            pcaVar(nFold, :, nTime)  = explained(1:numComps)/sum(explained);
            evTrialType(nFold, :, nTime) = mean(score(1:numTrials, :)).^2/totVar;
        end        
    end
        
    figure;
    subplot(2, 1, 1)
    hold on
    
    shadedErrorBar(DataSetList(nData).params.timeSeries, squeeze(mean(pcaVar(:, 1, :), 1)),...
            squeeze(std(pcaVar(:, 1, :), [], 1))/sqrt(numFold),...
            {'-', 'linewid', 1.0, 'color', 'k'}, 0.5);
        
    shadedErrorBar(DataSetList(nData).params.timeSeries, squeeze(mean(pcaVar(:, 2, :), 1)),...
            squeeze(std(pcaVar(:, 2, :), [], 1))/sqrt(numFold),...
            {'-', 'linewid', 1.0, 'color', 'm'}, 0.5);  
    legend({'PC1', 'PC2'})    
    axis xy;
    xlim([min(DataSetList(nData).params.timeSeries) max(DataSetList(nData).params.timeSeries)]);
    ylim([0 0.5])
    gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 0.5)
    box off;
    hold off;
    xlabel('Time (s)');
    ylabel('frac. PC var.');
    set(gca, 'TickDir', 'out')

    subplot(2, 1, 2)
    hold on
    shadedErrorBar(DataSetList(nData).params.timeSeries, squeeze(mean(evTrialType(:, 1, :), 1)),...
            squeeze(std(evTrialType(:, 1, :), [], 1))/sqrt(numFold),...
            {'-', 'linewid', 1.0, 'color', 'k'}, 0.5);
        
    shadedErrorBar(DataSetList(nData).params.timeSeries, squeeze(mean(evTrialType(:, 2, :), 1)),...
            squeeze(std(evTrialType(:, 2, :), [], 1))/sqrt(numFold),...
            {'-', 'linewid', 1.0, 'color', 'm'}, 0.5);  
    legend({'PC1', 'PC2'})    
    axis xy;
    xlim([min(DataSetList(nData).params.timeSeries) max(DataSetList(nData).params.timeSeries)]);
    ylim([0 0.5])
    gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 0.5)
    box off;
    hold off;
    xlabel('Time (s)');
    ylabel('frac. EV stim.');
    set(gca, 'TickDir', 'out')
    
    setPrint(2, 6*2, [PlotDir 'CollectedUnitsPCALDACorr/SimilarityPCA_500_' DataSetList(nData).name])
end

close all
