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
filterLength                  = 11;
filterStep                    = linspace(-filterLength / 2, filterLength / 2, filterLength);
filterInUse                   = exp(-filterStep .^ 2 / (2 * sigma ^ 2));
filterInUse                   = filterInUse / sum (filterInUse);

orgGroup                      = [4 4 5 5 6 6 7 7 7];
tarGroup                      = [1 2 1 3 2 3 1 2 3];


for nData      = 2:length(DataSetList)
    load([TempDatDir DataSetList(nData).name '.mat'])
    load([TempDatDir 'ModeledLogPValue_' DataSetList(nData).name '.mat'], 'logPValue', 'logPValueEpoch')
    unitGroup = getAnovaLogPSpikeEpoch (logPValue, logPValueEpoch, DataSetList(nData).params);
    
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
        
        for nNeuron     = 1:2
            
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
%             smoothedUnitData = nUnitData;
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
    
    setPrint(8*9, 6*6, [PlotDir 'ModeledComparingSingleUnitsAnova/SingleUnitsAnovaExampleNeuron_' DataSetList(nData).name], 'tif')
        
end

close all