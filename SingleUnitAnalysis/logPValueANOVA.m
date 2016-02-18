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

for nData      = 1:length(DataSetList)-1
%     load([TempDatDir DataSetList(nData).name '.mat'])
%     sigma                         = 0.10 / DataSetList(1).params.binsize; % 200 ms
%     filterLength                  = 10;
%     filterStep                    = linspace(-filterLength / 2, filterLength / 2, filterLength);
%     filterInUse                   = exp(-filterStep .^ 2 / (2 * sigma ^ 2));
%     filterInUse                   = filterInUse / sum (filterInUse);
% 
%     logPValue     = getLogPValueSpike(nDataSet, filterInUse);
%     save([TempDatDir 'LogPValue_' DataSetList(nData).name '.mat'], 'logPValue')
%     logPValueEpoch= getLogPValueSpikeEpoch(nDataSet, DataSetList(nData).params);
%     save([TempDatDir 'LogPValue_' DataSetList(nData).name '.mat'], 'logPValueEpoch', '-append')
    load([TempDatDir 'LogPValue_' DataSetList(nData).name '.mat'], 'logPValue', 'logPValueEpoch')
    unitGroup = plotAnovaLogPSpikeEpoch (logPValue, logPValueEpoch, DataSetList(nData).params);
    setPrint(10, 6*4, [PlotDir 'SingleUnitsAnova/SingleUnitsAnova_' DataSetList(nData).name], 'tiff')
    % plotAnovaLogPSpikeExampleNeurons (logPValue, DataSetList(nData), nDataSet, filterInUse);
    % % setPrint(20, 12*4, [PlotDir 'SingleUnitsAnovaExampleNeuron_' DataSetList(nData).name], 'tiff')
    
    depth                        = [DataSetList(nData).cellinfo(:).depth];    
    depthStart                   = 100;
    depthBin                     = 50;
    depthEnd                     = 900;    
    depth                        = floor((depth-depthStart)/depthBin)*depthBin+depthStart;
    depth(depth>depthEnd)        = depthEnd;
    depth(depth<depthStart)      = depthStart;
    uniqueDepth                  = depthStart:depthBin:depthEnd;
    depthStrings                 = cell(length(uniqueDepth),1);  
    depthStrings(1:end)          = {''};
    if length(uniqueDepth)       <=3
        depthStrings             = cellstr(num2str(uniqueDepth'));
    else
        stepLength               = floor(length(uniqueDepth)/3);
        depthStrings(1:stepLength:end) = cellstr(num2str(uniqueDepth(1:stepLength:end)'));
    end
    
    groupCounts = zeros(length(uniqueDepth), 7);
    
    for nGroup = 1:7
        nUnitGroup = unitGroup == nGroup;
        for nDepth = 1:length(uniqueDepth)
            groupCounts(nDepth, nGroup) = sum(nUnitGroup & depth' == uniqueDepth(nDepth));
        end
    end
    
%     zeroGroups  = sum(groupCounts, 2) == 0;
%     uniqueDepth = uniqueDepth(~zeroGroups);
%     groupCounts = groupCounts(~zeroGroups, :);
    groupPerCounts = bsxfun(@rdivide, groupCounts, sum(groupCounts, 2));
    
    figure;

    subplot(1, 2, 1)
    barh(uniqueDepth, groupPerCounts, 'stack', 'edgecolor', 'none');
    caxis([1 8])
    xlim([0 1])
    box off
    xlabel('% cell type')
    ylabel('Depth (um)')
    ylim([0 950])
    set(gca, 'yTick', 0:300:900)
    
    subplot(1, 2, 2)
    bar(uniqueDepth, sum(groupCounts,2))
    ylabel('# cells')
    xlabel('Depth (um)')
    xlim([0 950])
    set(gca, 'xTick', 0:300:900)
    
%     if length(uniqueDepth) > 1
%         barh(1:length(uniqueDepth), groupPerCounts, 'stack', 'edgecolor', 'none');
%         for nDepth = 1:length(uniqueDepth)
%             text(1.1, nDepth, ['Depth: ' num2str(uniqueDepth(nDepth)) '; # Cell: ' num2str(sum(groupCounts(nDepth,:)))],'fontsize',4);
%         end
%     else
%         barh(1:3, [zeros(1, 7); groupPerCounts; zeros(1, 7)], 'stack', 'edgecolor', 'none');
%         text(1.1, 2, ['Depth: ' num2str(uniqueDepth) '; # Cell: ' num2str(sum(groupCounts))],'fontsize',4);
%     end
%     colormap(parula(8))
%     caxis([1 8])
%     if length(uniqueDepth)>1
%         ylim([0.5 length(uniqueDepth)+0.5])
%     else
%         ylim([0.5 3.5])
%     end
%     xlim([0 2])
%     box off
%     axis off
%     ylabel('% cell type')
%     xlabel('Depth (um)')
    
    setPrint(8*2, 6, [PlotDir 'SingleUnitsAnova/SingleUnitsAnovaDepth_' DataSetList(nData).name], 'svg')
end

close all