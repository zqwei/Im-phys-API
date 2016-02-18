%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cell type categorization
% homogenous cell        : no change of selectivity
% heterogenous cell      : change of selectivity
% no-selective cell      : no selectivity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../Func');
setDir;
load ([TempDatDir 'DataListModeled.mat']);


if ~exist([PlotDir 'ModeledComparingSingleUnitsTscore'],'dir')
    mkdir([PlotDir 'ModeledComparingSingleUnitsTscore'])
end

nData          = 1;
load([TempDatDir DataSetList(nData).name '.mat'])
logPValueEpoch= getLogPValueTscoreSpikeEpoch(nDataSet, DataSetList(nData).params);
unitSpikeGroup = plotTtestLogPSpikeEpoch (logPValueEpoch);
spikeDataSet    = nDataSet;

sigma                         = 0.15 / DataSetList(1).params.binsize; % 200 ms
filterLength                  = 10;
filterStep                    = linspace(-filterLength / 2, filterLength / 2, filterLength);
filterInUse                   = exp(-filterStep .^ 2 / (2 * sigma ^ 2));
filterInUse                   = filterInUse / sum (filterInUse);


for nData      = 2:length(DataSetList)
    load([TempDatDir DataSetList(nData).name '.mat'])
    logPValueEpoch= getLogPValueTscoreSpikeEpoch(nDataSet, DataSetList(nData).params);
    unitGroup = plotTtestLogPSpikeEpoch (logPValueEpoch);    
%     groupCounts = zeros(3, 3);
%     
%     for nGroup  = 1:3
%         for mGroup = 1:3
%             groupCounts(nGroup, mGroup) = sum(unitGroup == mGroup-1 & unitSpikeGroup == nGroup-1);
%         end
%     end
%     
%     groupPerCounts = bsxfun(@rdivide, groupCounts, sum(groupCounts, 2));
%     
%     figure;
%     barh(1:3, groupPerCounts, 'stack', 'edgecolor', 'none');
%     caxis([1 3])
%     xlim([0 1])
%     box off
%     xlabel('% cell type')
%     ylabel('cell type')
%     ylim([0.5 3.5])
%     set(gca, 'yTick', 1:3, 'yTickLabel', {'Non.', 'Homo.', 'Dynamicial'})
%     setPrint(6, 4.5, [PlotDir 'ModeledComparingSingleUnitsTscore/SingleUnitsTscore_' DataSetList(nData).name], 'pdf')
    
    
    %%% Example neurons Dynamical --> Homo neuron
    figure;
    neuronIndice    = find(unitGroup == 1 & unitSpikeGroup == 2);
    
    color_index   = [0.7 0.0 0.0; % CL
                 0.0 0.0 0.7; % CR
                 1.0 0.7 0.7; % EL
                 0.7 0.7 1.0]; % ER

    
    for nExample = 1:10  
        neuronIndex    = neuronIndice(nExample);
        subplot(8, 5, mod(nExample-1, 5)+1 + floor((nExample-1)/5)*10)
        hold on;
        nUnitData        = spikeDataSet(neuronIndex).unit_yes_trial;
        smoothedUnitData = filter(filterInUse, 1, nUnitData, [], 2);
        shadedErrorBar(DataSetList(nData).params.timeSeries, mean(smoothedUnitData, 1),...
            std(smoothedUnitData, 1)/sqrt(size(smoothedUnitData, 1)),...
            {'-', 'linewid', 1.0, 'color', color_index(1, :)}, 0.5);
        nUnitData        = spikeDataSet(neuronIndex).unit_no_trial;
        smoothedUnitData = filter(filterInUse, 1, nUnitData, [], 2);
        shadedErrorBar(DataSetList(nData).params.timeSeries, mean(smoothedUnitData, 1),...
            std(smoothedUnitData, 1)/sqrt(size(smoothedUnitData, 1)),...
            {'-', 'linewid', 1.0, 'color', color_index(2, :)}, 0.5);
        gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 1.0)
        hold off;
        legend('boxoff');
        ylabel('Firing rate (Hz)');
        xlabel('Time (s)');
        xlim([DataSetList(nData).params.timeSeries(1) DataSetList(nData).params.timeSeries(end)]);
    
        subplot(8, 5, mod(nExample-1, 5)+1 + floor((nExample-1)/5)*10+5)
        hold on;
        nUnitData        = nDataSet(neuronIndex).unit_yes_trial;
        smoothedUnitData = filter(filterInUse, 1, nUnitData, [], 2);
        shadedErrorBar(DataSetList(nData).params.timeSeries, mean(smoothedUnitData, 1),...
            std(smoothedUnitData, 1)/sqrt(size(smoothedUnitData, 1)),...
            {'-', 'linewid', 1.0, 'color', color_index(1, :)}, 0.5);
        nUnitData        = nDataSet(neuronIndex).unit_no_trial;
        smoothedUnitData = filter(filterInUse, 1, nUnitData, [], 2);
        shadedErrorBar(DataSetList(nData).params.timeSeries, mean(smoothedUnitData, 1),...
            std(smoothedUnitData, 1)/sqrt(size(smoothedUnitData, 1)),...
            {'-', 'linewid', 1.0, 'color', color_index(2, :)}, 0.5);
        gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 1.0)
        hold off;
        legend('boxoff');
        ylabel('DF/F');
        xlabel('Time (s)');
        xlim([DataSetList(nData).params.timeSeries(1) DataSetList(nData).params.timeSeries(end)]);
    end
    
    %%% Example neurons Dynamical --> Dynamical neuron
        
    neuronIndice    = find(unitGroup == 2 & unitSpikeGroup == 2);
        
    for nExample = 1:10  
        neuronIndex    = neuronIndice(nExample);
        subplot(8, 5, mod(nExample-1, 5)+1 + floor((nExample-1)/5)*10 + 20)
        hold on;
        nUnitData        = spikeDataSet(neuronIndex).unit_yes_trial;
        smoothedUnitData = filter(filterInUse, 1, nUnitData, [], 2);
        shadedErrorBar(DataSetList(nData).params.timeSeries, mean(smoothedUnitData, 1),...
            std(smoothedUnitData, 1)/sqrt(size(smoothedUnitData, 1)),...
            {'-', 'linewid', 1.0, 'color', color_index(1, :)}, 0.5);
        nUnitData        = spikeDataSet(neuronIndex).unit_no_trial;
        smoothedUnitData = filter(filterInUse, 1, nUnitData, [], 2);
        shadedErrorBar(DataSetList(nData).params.timeSeries, mean(smoothedUnitData, 1),...
            std(smoothedUnitData, 1)/sqrt(size(smoothedUnitData, 1)),...
            {'-', 'linewid', 1.0, 'color', color_index(2, :)}, 0.5);
        gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 1.0)
        hold off;
        legend('boxoff');
        ylabel('Firing rate (Hz)');
        xlabel('Time (s)');
        xlim([DataSetList(nData).params.timeSeries(1) DataSetList(nData).params.timeSeries(end)]);
    
        subplot(8, 5, mod(nExample-1, 5)+1 + floor((nExample-1)/5)*10+ 25)
        hold on;
        nUnitData        = nDataSet(neuronIndex).unit_yes_trial;
        smoothedUnitData = filter(filterInUse, 1, nUnitData, [], 2);
        shadedErrorBar(DataSetList(nData).params.timeSeries, mean(smoothedUnitData, 1),...
            std(smoothedUnitData, 1)/sqrt(size(smoothedUnitData, 1)),...
            {'-', 'linewid', 1.0, 'color', color_index(1, :)}, 0.5);
        nUnitData        = nDataSet(neuronIndex).unit_no_trial;
        smoothedUnitData = filter(filterInUse, 1, nUnitData, [], 2);
        shadedErrorBar(DataSetList(nData).params.timeSeries, mean(smoothedUnitData, 1),...
            std(smoothedUnitData, 1)/sqrt(size(smoothedUnitData, 1)),...
            {'-', 'linewid', 1.0, 'color', color_index(2, :)}, 0.5);
        gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 1.0)
        hold off;
        legend('boxoff');
        ylabel('DF/F');
        xlabel('Time (s)');
        xlim([DataSetList(nData).params.timeSeries(1) DataSetList(nData).params.timeSeries(end)]);
    end
    
    setPrint(8*5, 6*8, [PlotDir 'ModeledComparingSingleUnitsTscore/SingleUnitsTscoreExampleNeuron_' DataSetList(nData).name], 'tif')
    
    
    
end

% close all