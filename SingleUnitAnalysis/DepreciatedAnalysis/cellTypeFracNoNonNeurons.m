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
    load([TempDatDir 'LogPValueTscore_' DataSetList(nData).name '.mat'], 'logPValue', 'logPValueEpoch')
    unitGroup = plotTtestLogPSpikeEpoch (logPValueEpoch);
    sizeGroup = histcounts(unitGroup, 0:3);
    figure('Visible', 'off');
    groupNames      = {'Non.', 'Homo.', 'Dynamical'};
    pie(sizeGroup(2:end))
    colormap(cmap(2:end, :))
    set(gca, 'TickDir', 'out')
    setPrint(8, 6, [PlotDir 'SingleUnitsTscore/SingleUnitsTscoreNoNon_' DataSetList(nData).name], 'svg')
end



if ~exist([PlotDir 'SingleUnitsAnova'],'dir')
    mkdir([PlotDir 'SingleUnitsAnova'])
end

cmap = cbrewer('qual', 'Set1', 10, 'cubic');

for nData      = 1:length(DataSetList)-1
    load([TempDatDir 'LogPValue_' DataSetList(nData).name '.mat'], 'logPValue', 'logPValueEpoch')
    unitGroup = getAnovaLogPSpikeEpoch(logPValueEpoch);
    plotAnovaLogPSpikeEpoch (logPValue, unitGroup, DataSetList(nData).params);    
    setPrint(10, 6*3, [PlotDir 'SingleUnitsAnova/SingleUnitsAnovaDynamic_' DataSetList(nData).name], 'svg')
    figure
    groupNames      = {'Pole', 'Lick', 'Reward', 'PL', 'PR', 'LR', 'PLR', 'Non.'};
    sizeGroup      = histcounts(unitGroup, 1:9);
    h = pie(sizeGroup(1:end-1));
    h = findobj(h, 'Type', 'patch');
    for nh = 1:length(h)
        set(h(nh), 'FaceColor', cmap(nh, :));
    end
    setPrint(8, 6, [PlotDir 'SingleUnitsAnova/SingleUnitsAnovaNoNon_' DataSetList(nData).name], 'svg')

end

close all