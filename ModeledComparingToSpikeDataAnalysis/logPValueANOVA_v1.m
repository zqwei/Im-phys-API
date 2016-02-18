%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cell type categorization
% AP cell                : sensory
% LR cell                : decision
% CE cell                : reward
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../Func');
setDir;
load ([TempDatDir 'DataListModeled.mat']);


if ~exist([PlotDir 'ModeledComparingSingleUnitsAnova'],'dir')
    mkdir([PlotDir 'ModeledComparingSingleUnitsAnova'])
end

nData         = 1;
load([TempDatDir 'ModeledLogPValue_' DataSetList(nData).name '.mat'], 'logPValue', 'logPValueEpoch')
unitSpikeGroup = plotAnovaLogPSpikeEpoch (logPValue, logPValueEpoch, DataSetList(nData).params);
close all

load([TempDatDir DataSetList(nData).name '.mat'])
spikeDataSet   = nDataSet;
minNumTrial    = 10;
numUnit        = length(spikeDataSet);


sigma                         = 0.15 / DataSetList(1).params.binsize; % 200 ms
filterLength                  = 10;
filterStep                    = linspace(-filterLength / 2, filterLength / 2, filterLength);
filterInUse                   = exp(-filterStep .^ 2 / (2 * sigma ^ 2));
filterInUse                   = filterInUse / sum (filterInUse);


for nData      = 2:length(DataSetList)
    load([TempDatDir DataSetList(nData).name '.mat'])
    load([TempDatDir 'ModeledLogPValue_' DataSetList(nData).name '.mat'], 'logPValue', 'logPValueEpoch')
    unitGroup = plotAnovaLogPSpikeEpoch (logPValue, logPValueEpoch, DataSetList(nData).params);
    close all
    
%     groupCounts = zeros(8, 8);
%     
%     for nGroup  = 1:8
%         for mGroup = 1:8
%             groupCounts(nGroup, mGroup) = sum(unitGroup == mGroup & unitSpikeGroup == nGroup);
%         end
%     end
%     
%     groupPerCounts = bsxfun(@rdivide, groupCounts, sum(groupCounts, 2));
%     
%     figure;
%     barh(1:8, groupPerCounts, 'stack', 'edgecolor', 'none');
%     caxis([1 8])
%     xlim([0 1])
%     box off
%     xlabel('% cell type')
%     ylabel('cell type')
%     ylim([0.5 8.5])
%     set(gca, 'yTick', 1:8, 'yTickLabel', {'P', 'L', 'R', 'PL', 'PR', 'LR', 'PLR', 'Non'})
%     setPrint(8, 6, [PlotDir 'ModeledComparingSingleUnitsAnova/SingleUnitsAnovaSwitch_' DataSetList(nData).name], 'pdf')
    
    %%% Example neurons PR --> P neuron
    figure;
    numTrialMat    = false(numUnit, 1);
    neuronIndex    = find(unitGroup == 1 & unitSpikeGroup == 5);
    
    for nUnit      = neuronIndex'
        numTrialMat(nUnit) =    length(nDataSet(nUnit).unit_yes_trial_index) > minNumTrial && ...
                        length(nDataSet(nUnit).unit_no_trial_index ) > minNumTrial && ...
                        length(nDataSet(nUnit).unit_yes_error_index) > minNumTrial && ...
                        length(nDataSet(nUnit).unit_no_error_index ) > minNumTrial;
    end
    
    neuronIndex    = find(numTrialMat);
    neuronIndex    = 168;%neuronIndex(2);
    
    subplot(3, 1, 1)
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
    
    subplot(3, 1, 2)
    hold on;
    nUnitData        = spikeDataSet(neuronIndex).unit_yes_trial;
    smoothedUnitData = filter(filterInUse, 1, nUnitData, [], 2);
    plot(DataSetList(nData).params.timeSeries, mean(smoothedUnitData, 1),'-', 'linewid', 1.0, 'color', color_index(1, :));
    nUnitData        = spikeDataSet(neuronIndex).unit_no_trial;
    smoothedUnitData = filter(filterInUse, 1, nUnitData, [], 2);
    plot(DataSetList(nData).params.timeSeries, mean(smoothedUnitData, 1),'-', 'linewid', 1.0, 'color', color_index(2, :));
    nUnitData        = spikeDataSet(neuronIndex).unit_yes_error;
    smoothedUnitData = filter(filterInUse, 1, nUnitData, [], 2);
    plot(DataSetList(nData).params.timeSeries, mean(smoothedUnitData, 1),'-', 'linewid', 1.0, 'color', color_index(3, :));
    nUnitData        = spikeDataSet(neuronIndex).unit_no_error;
    smoothedUnitData = filter(filterInUse, 1, nUnitData, [], 2);
    plot(DataSetList(nData).params.timeSeries, mean(smoothedUnitData, 1),'-', 'linewid', 1.0, 'color', color_index(4, :));
    gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 1.0)
    hold off;
    legend('CA', 'CP', 'EA', 'EP', 'Orientation', 'horizontal', 'Location', 'northoutside');
    legend('boxoff');
    ylabel('Firing rate (Hz)');
    xlabel('Time (s)');
    xlim([DataSetList(nData).params.timeSeries(1) DataSetList(nData).params.timeSeries(end)]);
    
    
    subplot(3, 1, 3)
    hold on;
    nUnitData        = nDataSet(neuronIndex).unit_yes_trial;
    smoothedUnitData = filter(filterInUse, 1, nUnitData, [], 2);
    plot(DataSetList(nData).params.timeSeries, mean(smoothedUnitData, 1),'-', 'linewid', 1.0, 'color', color_index(1, :));
    nUnitData        = nDataSet(neuronIndex).unit_no_trial;
    smoothedUnitData = filter(filterInUse, 1, nUnitData, [], 2);
    plot(DataSetList(nData).params.timeSeries, mean(smoothedUnitData, 1),'-', 'linewid', 1.0, 'color', color_index(2, :));
    nUnitData        = nDataSet(neuronIndex).unit_yes_error;
    smoothedUnitData = filter(filterInUse, 1, nUnitData, [], 2);
    plot(DataSetList(nData).params.timeSeries, mean(smoothedUnitData, 1),'-', 'linewid', 1.0, 'color', color_index(3, :));
    nUnitData        = nDataSet(neuronIndex).unit_no_error;
    smoothedUnitData = filter(filterInUse, 1, nUnitData, [], 2);
    plot(DataSetList(nData).params.timeSeries, mean(smoothedUnitData, 1),'-', 'linewid', 1.0, 'color', color_index(4, :));
    gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 1.0)
    hold off;
    legend('CA', 'CP', 'EA', 'EP', 'Orientation', 'horizontal', 'Location', 'northoutside');
    legend('boxoff');
    ylabel('DF/F');
    xlabel('Time (s)');
    xlim([DataSetList(nData).params.timeSeries(1) DataSetList(nData).params.timeSeries(end)]);
    
    setPrint(8, 6*3, [PlotDir 'ModeledComparingSingleUnitsAnova/SingleUnitsAnovaExampleNeuron_' DataSetList(nData).name], 'tif')
    
    %%% Example neurons LR --> L neuron
%     figure;
%     numTrialMat    = false(numUnit, 1);
%     neuronIndex    = find(unitGroup == 2 & unitSpikeGroup == 6);
%     
%     for nUnit      = neuronIndex'
%         numTrialMat(nUnit) =    length(nDataSet(nUnit).unit_yes_trial_index) > minNumTrial && ...
%                         length(nDataSet(nUnit).unit_no_trial_index ) > minNumTrial && ...
%                         length(nDataSet(nUnit).unit_yes_error_index) > minNumTrial && ...
%                         length(nDataSet(nUnit).unit_no_error_index ) > minNumTrial;
%     end
    
%     neuronIndex    = find(numTrialMat, 2);
%     
%     subplot(3, 1, 1)
%     spkTimes{1}    = spikeDataSet(neuronIndex).unit_yes_trial_spk_time;
%     spkTimes{2}    = spikeDataSet(neuronIndex).unit_no_trial_spk_time;
%     spkTimes{3}    = spikeDataSet(neuronIndex).unit_yes_error_spk_time;
%     spkTimes{4}    = spikeDataSet(neuronIndex).unit_no_error_spk_time;    
%     
%     color_index   = [0.7 0.0 0.0; % CL
%                  0.0 0.0 0.7; % CR
%                  0.3 0.0 0.0; % EL
%                  0.0 0.0 0.3]; % ER
%     
%     for nPlot            = 1:4
%         hold on;
%         for nTrial       = 1:minNumTrial
%             spikeTimes   = spkTimes{nPlot}{nTrial};
%             spikeTrial   = ones(length(spikeTimes), 1) * (nTrial + (nPlot-1)*(minNumTrial+3));
%             plot(spikeTimes, spikeTrial, '.', 'color', color_index(nPlot, :));
%         end
%     end
%     
%     gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 1.0);
%     hold off;
%     ylim([1 (minNumTrial + 3) * 4])
%     xlim([DataSetList(nData).params.timeSeries(1) DataSetList(nData).params.timeSeries(end)]);
%     axis off
%     
%     subplot(3, 1, 2)
%     hold on;
%     nUnitData        = spikeDataSet(neuronIndex).unit_yes_trial;
%     smoothedUnitData = filter(filterInUse, 1, nUnitData, [], 2);
%     plot(DataSetList(nData).params.timeSeries, mean(smoothedUnitData, 1),'-', 'linewid', 1.0, 'color', color_index(1, :));
%     nUnitData        = spikeDataSet(neuronIndex).unit_no_trial;
%     smoothedUnitData = filter(filterInUse, 1, nUnitData, [], 2);
%     plot(DataSetList(nData).params.timeSeries, mean(smoothedUnitData, 1),'-', 'linewid', 1.0, 'color', color_index(2, :));
%     nUnitData        = spikeDataSet(neuronIndex).unit_yes_error;
%     smoothedUnitData = filter(filterInUse, 1, nUnitData, [], 2);
%     plot(DataSetList(nData).params.timeSeries, mean(smoothedUnitData, 1),'-', 'linewid', 1.0, 'color', color_index(3, :));
%     nUnitData        = spikeDataSet(neuronIndex).unit_no_error;
%     smoothedUnitData = filter(filterInUse, 1, nUnitData, [], 2);
%     plot(DataSetList(nData).params.timeSeries, mean(smoothedUnitData, 1),'-', 'linewid', 1.0, 'color', color_index(4, :));
%     gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 1.0)
%     hold off;
%     legend('CA', 'CP', 'EA', 'EP', 'Orientation', 'horizontal', 'Location', 'northoutside');
%     legend('boxoff');
%     ylabel('Firing rate (Hz)');
%     xlabel('Time (s)');
%     xlim([DataSetList(nData).params.timeSeries(1) DataSetList(nData).params.timeSeries(end)]);
%     
%     
%     subplot(3, 1, 3)
%     hold on;
%     nUnitData        = nDataSet(neuronIndex).unit_yes_trial;
%     smoothedUnitData = filter(filterInUse, 1, nUnitData, [], 2);
%     plot(DataSetList(nData).params.timeSeries, mean(smoothedUnitData, 1),'-', 'linewid', 1.0, 'color', color_index(1, :));
%     nUnitData        = nDataSet(neuronIndex).unit_no_trial;
%     smoothedUnitData = filter(filterInUse, 1, nUnitData, [], 2);
%     plot(DataSetList(nData).params.timeSeries, mean(smoothedUnitData, 1),'-', 'linewid', 1.0, 'color', color_index(2, :));
%     nUnitData        = nDataSet(neuronIndex).unit_yes_error;
%     smoothedUnitData = filter(filterInUse, 1, nUnitData, [], 2);
%     plot(DataSetList(nData).params.timeSeries, mean(smoothedUnitData, 1),'-', 'linewid', 1.0, 'color', color_index(3, :));
%     nUnitData        = nDataSet(neuronIndex).unit_no_error;
%     smoothedUnitData = filter(filterInUse, 1, nUnitData, [], 2);
%     plot(DataSetList(nData).params.timeSeries, mean(smoothedUnitData, 1),'-', 'linewid', 1.0, 'color', color_index(4, :));
%     gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 1.0)
%     hold off;
%     legend('CA', 'CP', 'EA', 'EP', 'Orientation', 'horizontal', 'Location', 'northoutside');
%     legend('boxoff');
%     ylabel('DF/F');
%     xlabel('Time (s)');
%     xlim([DataSetList(nData).params.timeSeries(1) DataSetList(nData).params.timeSeries(end)]);

    
end

% close all