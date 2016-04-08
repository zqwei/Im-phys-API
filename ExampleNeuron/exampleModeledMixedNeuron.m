%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cell type categorization
% AP cell                : sensory
% LR cell                : decision
% CE cell                : reward
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../Func');
setDir;
load ([TempDatDir 'DataListShuffle.mat']);
nData         = 1;
load([TempDatDir 'LogPValue_' DataSetList(nData).name '.mat'], 'logPValueEpoch')
unitSpikeGroup = getAnovaLogPSpikeEpoch(logPValueEpoch);
load([TempDatDir DataSetList(nData).name '.mat'])
spikeDataSet   = nDataSet;
minNumTrial    = 10;
numUnit        = length(spikeDataSet);
sigma                         = 0.15 / DataSetList(1).params.binsize; % 200 ms
filterLength                  = 11;
filterStep                    = linspace(-filterLength / 2, filterLength / 2, filterLength);
filterInUse                   = exp(-filterStep .^ 2 / (2 * sigma ^ 2));
filterInUse                   = filterInUse / sum (filterInUse);


load ([TempDatDir 'DataListS2CModel.mat']);
if ~exist([PlotDir 'S2CModel'],'dir')
    mkdir([PlotDir 'S2CModel'])
end
cmap = cbrewer('qual', 'Set1', 10, 'cubic');


orgGroup                      = [4 4 5 5 6 6 7 7 7];
tarGroup                      = [1 2 1 3 2 3 1 2 3];


for nData      = 1:length(DataSetList)
    load([TempDatDir DataSetList(nData).name '.mat'])
    
    logPValueEpoch= getLogPValueSpikeEpoch(nDataSet, DataSetList(nData).params);
    unitGroup = getAnovaLogPSpikeEpoch(logPValueEpoch);
    
    figure;
    groupCounts = zeros(8, 8);
    
    for nGroup  = 1:8
        for mGroup = 1:8
            groupCounts(nGroup, mGroup) = sum(unitGroup == mGroup & unitSpikeGroup == nGroup);
        end
    end
    
    groupPerCounts = bsxfun(@rdivide, groupCounts, sum(groupCounts, 2));
    barh(1:8, groupPerCounts, 'stack', 'edgecolor', 'none');
    caxis([1 8])
    xlim([0 1])
    box off
    xlabel('frac. cell type changed in model ')
    ylabel('original cell type')
    ylim([0.5 8.5])
    colormap(cmap(1:8, :))
    set(gca, 'TickDir', 'out')
    set(gca, 'yTick', 1:8, 'yTickLabel', {'P', 'L', 'R', 'PL', 'PR', 'LR', 'PLR', 'Non'})
    setPrint(8, 6, [PlotDir 'S2CModel/SingleUnitsAnovaTrans_' DataSetList(nData).name])
    
    figure;
    for nExample   = 1:length(orgGroup)

        numTrialMat    = false(numUnit, 1);
        neuronIndex    = find(unitGroup == tarGroup(nExample) & unitSpikeGroup == orgGroup(nExample));

        for nUnit      = neuronIndex'
            numTrialMat(nUnit) =    length(nDataSet(nUnit).unit_yes_trial_index) > minNumTrial && ...
                            length(nDataSet(nUnit).unit_no_trial_index ) > minNumTrial && ...
                            length(nDataSet(nUnit).unit_yes_error_index) > minNumTrial && ...
                            length(nDataSet(nUnit).unit_no_error_index ) > minNumTrial;
        end

        neuronIndice    = find(numTrialMat);
        
        numNeuron       = min(2, length(neuronIndice));
        
        for nNeuron     = 1:numNeuron
            
            neuronIndex = neuronIndice(nNeuron);
            
            subplot(6, 9, nExample + ((nNeuron-1)*3)*9 )
            spkTimes{1}    = spikeDataSet(neuronIndex).unit_yes_trial_spk_time;
            spkTimes{2}    = spikeDataSet(neuronIndex).unit_no_trial_spk_time;
            spkTimes{3}    = spikeDataSet(neuronIndex).unit_yes_error_spk_time;
            spkTimes{4}    = spikeDataSet(neuronIndex).unit_no_error_spk_time;    

            color_index   = [0.7 0.0 0.0; % CL
                         0.0 0.0 0.7; % CR
                         1.0 0.7 0.7; % EL
                         0.7 0.7 1.0]; % ER

            spkLoc         = 0;
            for nPlot            = 1:4
                hold on;
                for nTrial       = 1:length(spkTimes{nPlot})
                    spikeTimes   = spkTimes{nPlot}{nTrial};
                    spikeTrial   = ones(length(spikeTimes), 1) * (nTrial + spkLoc);
                    plot(spikeTimes, spikeTrial, '.', 'color', color_index(nPlot, :));
                end
                spkLoc     = spkLoc + length(spkTimes{nPlot}) + 3;
            end
            
            

            gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 1.0);
            hold off;
            ylim([1 spkLoc])
            xlim([DataSetList(nData).params.timeSeries(1) DataSetList(nData).params.timeSeries(end)]);
            axis off

            subplot(6, 9, nExample + ((nNeuron-1)*3+1)*9)
            hold on;
            nUnitData        = spikeDataSet(neuronIndex).unit_yes_trial;
            smoothedUnitData = getGaussianPSTH (filterInUse, nUnitData, 2);
            shadedErrorBar(DataSetList(nData).params.timeSeries, mean(smoothedUnitData, 1),...
            std(smoothedUnitData, 1)/sqrt(size(smoothedUnitData, 1)),...
            {'k-', 'linewid', 1.0, 'color', color_index(1, :)}, 0.5);
            nUnitData        = spikeDataSet(neuronIndex).unit_no_trial;
            smoothedUnitData = getGaussianPSTH (filterInUse, nUnitData, 2);
            shadedErrorBar(DataSetList(nData).params.timeSeries, mean(smoothedUnitData, 1),...
            std(smoothedUnitData, 1)/sqrt(size(smoothedUnitData, 1)),...
            {'-', 'linewid', 1.0, 'color', color_index(2, :)}, 0.5);
            nUnitData        = spikeDataSet(neuronIndex).unit_yes_error;
            smoothedUnitData = getGaussianPSTH (filterInUse, nUnitData, 2);
            shadedErrorBar(DataSetList(nData).params.timeSeries, mean(smoothedUnitData, 1),...
            std(smoothedUnitData, 1)/sqrt(size(smoothedUnitData, 1)),...
            {'-', 'linewid', 1.0, 'color', color_index(3, :)}, 0.5);
            nUnitData        = spikeDataSet(neuronIndex).unit_no_error;
            smoothedUnitData = getGaussianPSTH (filterInUse, nUnitData, 2);
            shadedErrorBar(DataSetList(nData).params.timeSeries, mean(smoothedUnitData, 1),...
            std(smoothedUnitData, 1)/sqrt(size(smoothedUnitData, 1)),...
            {'-', 'linewid', 1.0, 'color', color_index(4, :)}, 0.5);
            gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 1.0)
            hold off;
            ylabel('Firing rate (Hz)');
            xlabel('Time (s)');
            xlim([DataSetList(nData).params.timeSeries(1) DataSetList(nData).params.timeSeries(end)]);


            subplot(6, 9, nExample + ((nNeuron-1)*3+2)*9)
            hold on;
            nUnitData        = nDataSet(neuronIndex).unit_yes_trial;
            smoothedUnitData = getGaussianPSTH (filterInUse, nUnitData, 2);
            shadedErrorBar(DataSetList(nData).params.timeSeries, mean(smoothedUnitData, 1),...
            std(smoothedUnitData, 1)/sqrt(size(smoothedUnitData, 1)),...
            {'-', 'linewid', 1.0, 'color', color_index(1, :)}, 0.5);
            nUnitData        = nDataSet(neuronIndex).unit_no_trial;
            smoothedUnitData = getGaussianPSTH (filterInUse, nUnitData, 2);
            shadedErrorBar(DataSetList(nData).params.timeSeries, mean(smoothedUnitData, 1),...
            std(smoothedUnitData, 1)/sqrt(size(smoothedUnitData, 1)),...
            {'-', 'linewid', 1.0, 'color', color_index(2, :)}, 0.5);
            nUnitData        = nDataSet(neuronIndex).unit_yes_error;
            smoothedUnitData = getGaussianPSTH (filterInUse, nUnitData, 2);
            shadedErrorBar(DataSetList(nData).params.timeSeries, mean(smoothedUnitData, 1),...
            std(smoothedUnitData, 1)/sqrt(size(smoothedUnitData, 1)),...
            {'-', 'linewid', 1.0, 'color', color_index(3, :)}, 0.5);
            nUnitData        = nDataSet(neuronIndex).unit_no_error;
            smoothedUnitData = getGaussianPSTH (filterInUse, nUnitData, 2);
            shadedErrorBar(DataSetList(nData).params.timeSeries, mean(smoothedUnitData, 1),...
            std(smoothedUnitData, 1)/sqrt(size(smoothedUnitData, 1)),...
            {'-', 'linewid', 1.0, 'color', color_index(4, :)}, 0.5);
            gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 1.0)
            hold off;
            ylabel('DF/F');
            xlabel('Time (s)');
            xlim([DataSetList(nData).params.timeSeries(1) DataSetList(nData).params.timeSeries(end)]);
            
        end
        
    end
    
    setPrint(8*9, 6*6, [PlotDir 'S2CModel/SingleUnitsAnovaExampleNeuron_' DataSetList(nData).name])
        
end

close all