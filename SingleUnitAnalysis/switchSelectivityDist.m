%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cell type categorization
% homogenous cell        : no change of selectivity
% heterogenous cell      : change of selectivity
% no-selective cell      : no selectivity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../Func');
setDir;
load ([TempDatDir 'DataListShuffle.mat']);


if ~exist([PlotDir 'SingleUnitsTscore'],'dir')
    mkdir([PlotDir 'SingleUnitsTscore'])
end

cmap = cbrewer('qual', 'Set1', 3, 'cubic');

for nData      = 1:length(DataSetList)-1
%     load([TempDatDir DataSetList(nData).name '.mat'])
%     logPValueEpoch= getLogPValueTscoreSpikeEpoch(nDataSet, DataSetList(nData).params);
%     save([TempDatDir 'LogPValueTscore_' DataSetList(nData).name '.mat'], 'logPValueEpoch')
    load([TempDatDir 'LogPValueTscore_' DataSetList(nData).name '.mat'], 'logPValue', 'logPValueEpoch')
    unitGroup = plotTtestLogPSpikeEpoch (logPValueEpoch);
    sizeGroup = histcounts(unitGroup, 0:3);
    figure('Visible', 'off');
    groupNames      = {'Non.', 'Homo.', 'Dynamical'};
%     pie(sizeGroup, groupNames)
    pie(sizeGroup)
%     title('Distribution of cell types')
    colormap(cmap)
    set(gca, 'TickDir', 'out')
    setPrint(8, 6, [PlotDir 'SingleUnitsTscore/SingleUnitsTscore_' DataSetList(nData).name], 'svg')
    
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
    
    groupCounts = zeros(length(uniqueDepth), 3);
    
    for nGroup = 1:3
        nUnitGroup = unitGroup == nGroup-1;
        for nDepth = 1:length(uniqueDepth)
            groupCounts(nDepth, nGroup) = sum(nUnitGroup & depth' == uniqueDepth(nDepth));
        end
    end
    
    groupPerCounts = bsxfun(@rdivide, groupCounts, sum(groupCounts, 2));
    
    figure('Visible', 'off');
    subplot(1, 2, 1)
    barh(uniqueDepth, groupPerCounts, 'stack', 'edgecolor', 'none');
%     caxis([0 3])
    xlim([0 1])
    box off
    xlabel('% cell type')
    ylabel('Depth (um)')
    ylim([0 950])
    colormap(cmap)
    set(gca, 'yTick', 0:300:900)
    
    subplot(1, 2, 2)
    barh(uniqueDepth, sum(groupCounts,2),'k')
    xlabel('# cells')
    ylabel('Depth (um)')
    ylim([0 950])
    set(gca, 'yTick', 0:300:900)
    set(gca, 'TickDir', 'out')
    setPrint(8*2, 6, [PlotDir 'SingleUnitsTscore/SingleUnitsTscoreDepth_' DataSetList(nData).name], 'svg')
end

figure;
hold on
for nColor = 1:length(groupNames)
    plot(0, nColor, 's', 'color', cmap(nColor,:), 'MarkerFaceColor',cmap(nColor,:),'MarkerSize', 8)
    text(1, nColor, groupNames{nColor})
end
xlim([0 10])
hold off
axis off
setPrint(3, 2, [PlotDir 'SingleUnitsTscore/SingleUnitsTscore_Label'])
close all