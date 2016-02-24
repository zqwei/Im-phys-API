%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cell type categorization
% AP cell                : sensory
% LR cell                : decision
% CE cell                : reward
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../Func');
setDir;
load ([TempDatDir 'DataListShuffle.mat']);


if ~exist([PlotDir 'SingleUnitsAnova'],'dir')
    mkdir([PlotDir 'SingleUnitsAnova'])
end

numExample             = 1;

color_index   = [0.7 0.0 0.0; % CL
                 0.0 0.0 0.7; % CR
                 1.0 0.7 0.7; % EL
                 0.7 0.7 1.0]; % ER

nData         = 1;
load([TempDatDir DataSetList(nData).name '.mat'])             
sigma                         = 0.15 / DataSetList(1).params.binsize; % 200 ms
filterLength                  = 11;
filterStep                    = linspace(-filterLength / 2, filterLength / 2, filterLength);
filterInUse                   = exp(-filterStep .^ 2 / (2 * sigma ^ 2));
filterInUse                   = filterInUse / sum (filterInUse);  
textYLabels                   = {'Firing rate (Hz)', 'DF/F', 'DF/F', 'DF/F'};
groupNames                    = {'Pole', 'Lick', 'Reward', 'PL', 'PR', 'LR', 'PLR', 'Non.'};
minUnits                      = [15 2 3 15];

for nData              = [1 3 4]
    load([TempDatDir 'LogPValue_' DataSetList(nData).name '.mat'], 'logPValueEpoch')
    load([TempDatDir DataSetList(nData).name '.mat'])             

    unitGroup = getAnovaLogPSpikeEpoch(logPValueEpoch);

    figure;
    
    for nGroup         = 1:7
        neuronIndice   = find(unitGroup==nGroup);
        nExample       = 1;
        neuronIndex    = neuronIndice(nExample);
        maxExample     = length(neuronIndice);
        while size(nDataSet(neuronIndex).unit_no_error,1)<minUnits(nData) || size(nDataSet(neuronIndex).unit_yes_error,1)<minUnits(nData)
            if nExample < maxExample
                nExample    = nExample + 1;
                neuronIndex = neuronIndice(nExample);
            end
        end
        
        subplot(2, 4, nGroup)
        hold on;
        nUnitData        = nDataSet(neuronIndex).unit_yes_trial;
        smoothedUnitData = getGaussianPSTH (filterInUse, nUnitData, 2);
        shadedErrorBar(DataSetList(nData).params.timeSeries, mean(smoothedUnitData, 1),...
        std(smoothedUnitData, 1)/sqrt(size(smoothedUnitData, 1)),...
        {'k-', 'linewid', 1.0, 'color', color_index(1, :)}, 0.5);
        
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
        ylabel(textYLabels{nData});
        title(groupNames{nGroup})
        xlabel('Time (s)');
        xlim([DataSetList(nData).params.timeSeries(1) DataSetList(nData).params.timeSeries(end)]);
        set(gca, 'TickDir', 'out')
    end
    setPrint(8*4, 6*2, [PlotDir 'SingleUnitsAnova/SingleUnitsAnovaExampleNeuron_' DataSetList(nData).name])
end