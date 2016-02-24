%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Collected population PCA variance and EV of trial type over time
%
% fixed number of neurons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../Func');
setDir;
load ([TempDatDir 'DataListShuffle.mat']);

numUnits       = 400;
numTrials      = numUnits * 3;
numComps       = 15;
ROCThres       = 0.5;
trialType      = [true(numTrials, 1); false(numTrials, 1)];

if ~exist([PlotDir 'CollectedUnitsPCA'],'dir')
    mkdir([PlotDir 'CollectedUnitsPCA'])
end

for nData              = 1:length(DataSetList)-1
    load([TempDatDir DataSetList(nData).name '.mat']);
    numT               = size(nDataSet(1).unit_yes_trial, 2);
    pcaVar             = nan(numComps, numT);
    evTrialType        = nan(numComps, numT);

    selectedNeuronalIndex = DataSetList(nData).ActiveNeuronIndex';
    selectedNeuronalIndex = selectedHighROCneurons(nDataSet, DataSetList(nData).params, ROCThres, selectedNeuronalIndex);
    nDataSet              = nDataSet(selectedNeuronalIndex);
    
    if length(nDataSet)>numUnits    
        
        randPickUnits         = randperm(length(nDataSet));
        randPickUnits         = randPickUnits(1:numUnits);
        nDataSet              = nDataSet(randPickUnits);

        firingRates        = generatePCAData(nDataSet, numTrials);


        for nTime          = 1:numT
            nFiringRates       = squeeze(firingRates(nTime, :, :))';
            [~,score,latent,~,explained,~] = pca(nFiringRates, 'NumComponents', numComps);
            totVar            = sum(latent);
            pcaVar(:, nTime)  = explained(1:numComps)/sum(explained);
            evTrialType(:, nTime) = mean(score(1:numTrials, :)).^2/totVar;
        end

        save([TempDatDir 'PCATimeFixed' num2str(numUnits) 'Units_' DataSetList(nData).name '.mat'], 'pcaVar', 'evTrialType')

        figure;
        subplot(1, 2, 1)
        hold on
        imagesc(DataSetList(nData).params.timeSeries, 1:numComps, pcaVar)
        axis xy;
        xlim([min(DataSetList(nData).params.timeSeries) max(DataSetList(nData).params.timeSeries)]);
        ylim([1 numComps])
        gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 0.5)
        box off;
        hold off;
        xlabel('Time (s)');
        ylabel('Component index');
        title('PC Variance')
        colorbar

        subplot(1, 2, 2)
        hold on
        imagesc(DataSetList(nData).params.timeSeries, 1:numComps, evTrialType)
        axis xy;
        xlim([min(DataSetList(nData).params.timeSeries) max(DataSetList(nData).params.timeSeries)]);
        ylim([1 numComps])
        gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 0.5)
        box off;
        hold off;
        xlabel('Time (s)');
        ylabel('Component index');
        title('PC Trial Type EV')
        colorbar

        setPrint(8*2, 6, [PlotDir 'CollectedUnitsPCA/CollectedUnitsPCAFixed' num2str(numUnits) 'Units_' DataSetList(nData).name], 'pdf')
    end
end

close all
