%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cell type categorization
% homogenous cell        : no change of selectivity
% heterogenous cell      : change of selectivity
% no-selective cell      : no selectivity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../Func');
setDir;
load ([TempDatDir 'DataListModeled.mat']);


if ~exist([PlotDir 'ModeledSingleUnitsTscore'],'dir')
    mkdir([PlotDir 'ModeledSingleUnitsTscore'])
end

for nData      = 1:length(DataSetList)
    load([TempDatDir DataSetList(nData).name '.mat'])
    logPValueEpoch= getLogPValueTscoreSpikeEpoch(nDataSet, DataSetList(nData).params);
    save([TempDatDir 'LogPValueTscore_' DataSetList(nData).name '.mat'], 'logPValueEpoch')
    load([TempDatDir 'LogPValueTscore_' DataSetList(nData).name '.mat'], 'logPValue', 'logPValueEpoch')
    unitGroup = plotTtestLogPSpikeEpoch (logPValueEpoch);
    sizeGroup = histcounts(unitGroup, 0:3);
    figure;
    groupNames      = {'Non.', 'Homo.', 'Dynamicial'};
    pie(sizeGroup, groupNames)
    title('Distribution of cell types')
    setPrint(6, 4.5, [PlotDir 'ModeledSingleUnitsTscore/SingleUnitsTscore_' DataSetList(nData).name], 'pdf')
    
    depth                        = [DataSetList(nData).cellinfo(:).depth];    
    depthStart                   = 100;
    depthBin                     = 50;
    depthEnd                     = 900;    
    depth                        = floor((depth-depthStart)/depthBin)*depthBin+depthStart;
    depth(depth>depthEnd)        = depthEnd;
    uniqueDepth                  = depthStart:depthBin:depthEnd;
    depthStrings                 = cell(length(uniqueDepth),1);  
    depthStrings(1:end)          = {''};
    if length(uniqueDepth)       <=3
        depthStrings             = cellstr(num2str(uniqueDepth'));
    else
        stepLength               = floor(length(uniqueDepth)/3);
        depthStrings(1:stepLength:end) = cellstr(num2str(uniqueDepth(1:stepLength:end)'));
    end
    
    groupCounts = zeros(length(uniqueDepth), 3);
    
    for nGroup = 1:3
        nUnitGroup = unitGroup == nGroup-1;
        for nDepth = 1:length(uniqueDepth)
            groupCounts(nDepth, nGroup) = sum(nUnitGroup & depth' == uniqueDepth(nDepth));
        end
    end
    
    groupPerCounts = bsxfun(@rdivide, groupCounts, sum(groupCounts, 2));
    
    figure;
    subplot(1, 2, 1)
    barh(uniqueDepth, groupPerCounts, 'stack', 'edgecolor', 'none');
%     caxis([0 3])
    xlim([0 1])
    box off
    xlabel('% cell type')
    ylabel('Depth (um)')
    ylim([0 950])
    set(gca, 'yTick', 0:300:900)
    
    subplot(1, 2, 2)
    barh(uniqueDepth, sum(groupCounts,2))
    xlabel('# cells')
    ylabel('Depth (um)')
    ylim([0 950])
    set(gca, 'yTick', 0:300:900)
    
    setPrint(8*2, 6, [PlotDir 'ModeledSingleUnitsTscore/SingleUnitsTscoreDepth_' DataSetList(nData).name], 'pdf')
end

close all