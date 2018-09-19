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

cmap = cbrewer('qual', 'Set1', 9, 'cubic');
cmap = cmap([3, 5, 9], :);
groupColors = {cmap(3, :), cmap(2, :), cmap(1, :)};


for nData      = [1 10]
    if ~exist([TempDatDir DataSetList(nData).name '_withOLRemoval.mat'], 'file')
        load([TempDatDir DataSetList(nData).name '.mat']);
        neuronRemoveList = false(length(nDataSet), 1);
    else
        load([TempDatDir DataSetList(nData).name '_withOLRemoval.mat']);
    end
    
    depth_list          = [nDataSet.depth_in_um]';    
    oldDataSet          = nDataSet;
    selectedNeuronalIndex = depth_list < 471;
    
    nDataSet              = oldDataSet(selectedNeuronalIndex);
    unitGroup = getLogPValueTscoreSpikeTime(nDataSet, DataSetList(nData).params);
    sizeGroup = histcounts(unitGroup, 0:3);
    disp(sizeGroup/length(nDataSet))
    disp(length(nDataSet))
    figure('Visible', 'off');
    groupNames      = {['Non' newline 'n = ' num2str(sum(unitGroup==0))], ... 
                       ['Mono' newline 'n = ' num2str(sum(unitGroup==1))], ...
                       ['Multi' newline 'n = ' num2str(sum(unitGroup==2))]};
%     donut(sizeGroup, groupNames, groupColors);
    pie(sizeGroup, groupNames)
    axis off
%     legend('Location','eastoutside')
%     legend('boxoff')
    setPrint(8, 6, [PlotDir 'SingleUnitsTscore/SingleUnitsTscoreTime_' DataSetList(nData).name])
%     setPrint(8, 6, [PlotDir 'SingleUnitsTscore/' DataSetList(nData).name '_selectivity'], 'svg')

    depth                        = [nDataSet.depth_in_um]';
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
            groupCounts(nDepth, nGroup) = sum(nUnitGroup & depth == uniqueDepth(nDepth));
        end
    end

    [tab, chi2, p] = crosstab(unitGroup(unitGroup>0), depth(unitGroup>0)');
    disp(p)

    groupPerCounts = bsxfun(@rdivide, groupCounts, sum(groupCounts, 2));

    figure('Visible', 'off');
    subplot(1, 2, 1)
    barh(-uniqueDepth, groupPerCounts, 'stack', 'edgecolor', 'none');
    xlim([0 1])
    box off
    xlabel('% cell type')
    ylabel('Depth (um)')
    set(gca, 'TickDir', 'out')
    ylim([-850 0])
    set(gca, 'yTick', -800:400:0)
    colormap(cmap)

    subplot(1, 2, 2)
    barh(-uniqueDepth, sum(groupCounts,2),'k')
    xlabel('# cells')
    ylabel('Depth (um)')
    ylim([-850 0])
    set(gca, 'yTick', -800:400:0)
    set(gca, 'TickDir', 'out')
    setPrint(8*2, 6, [PlotDir 'SingleUnitsTscore/SingleUnitsTscoreTimeDepth_' DataSetList(nData).name])
end