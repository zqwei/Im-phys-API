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

% cmap = cbrewer('qual', 'Set1', 3, 'cubic');

cmap = [0.8000    0.8000    0.8000;
       1.0000    0.6000         0;
       0    0.8000         0];

for nData      = [3 4]
    load([TempDatDir DataSetList(nData).name '_withOLRemoval.mat'])
    unitGroup = getLogPValueTscoreSpikeTime(nDataSet, DataSetList(nData).params);    
    sizeGroup = histcounts(unitGroup, 0:3);
%     sizeGroup
%     sum(sizeGroup)
%     disp(sizeGroup(2)/sizeGroup(3))
    figure('Visible', 'off');
    groupNames      = {'Non.', 'Homo.', 'Dynamical'};
    pie(sizeGroup)
    colormap(cmap)
    set(gca, 'TickDir', 'out')
    setPrint(8, 6, [PlotDir 'SingleUnitsTscore/SingleUnitsTscoreTime_' DataSetList(nData).name '_withOLRemoval'])
    
    depth                        = [DataSetList(nData).cellinfo(:).depth];
    depth                        = depth(~neuronRemoveList);
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
    
%     [tab, chi2, p] = crosstab(unitGroup(unitGroup>0), depth(unitGroup>0)')
    
    groupPerCounts = bsxfun(@rdivide, groupCounts, sum(groupCounts, 2));
    
    figure('Visible', 'off');
    subplot(1, 2, 1)
    barh(-uniqueDepth, groupPerCounts, 'stack', 'edgecolor', 'none');
    xlim([0 1])
    box off
    xlabel('% cell type')
    ylabel('Depth (um)')
    set(gca, 'TickDir', 'out')
    ylim([-950 0])
    colormap(cmap)
    set(gca, 'yTick', -900:300:0)
    
    subplot(1, 2, 2)
    barh(-uniqueDepth, sum(groupCounts,2),'k')
    xlabel('# cells')
    ylabel('Depth (um)')
    ylim([-950 0])
    set(gca, 'yTick', -900:300:0)
    set(gca, 'TickDir', 'out')
    setPrint(8*2, 6, [PlotDir 'SingleUnitsTscore/SingleUnitsTscoreTimeDepth_' DataSetList(nData).name '_withOLRemoval'])
end