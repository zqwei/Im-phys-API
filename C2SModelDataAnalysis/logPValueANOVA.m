%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cell type categorization
% AP cell                : sensory
% LR cell                : decision
% CE cell                : reward
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../Func');
setDir;
load ([TempDatDir 'DataListC2SModel.mat']);


if ~exist([PlotDir 'C2SModel'],'dir')
    mkdir([PlotDir 'C2SModel'])
end

cmap = cbrewer('qual', 'Set1', 10, 'cubic');

for nData      = 1:length(DataSetList)
    load([TempDatDir DataSetList(nData).name '.mat'])
    logPValueEpoch= getLogPValueSpikeEpoch(nDataSet, DataSetList(nData).params);
    unitGroup = getAnovaLogPSpikeEpoch(logPValueEpoch);
    figure
    groupNames      = {'Pole', 'Lick', 'Reward', 'PL', 'PR', 'LR', 'PLR', 'Non.'};
    sizeGroup      = histcounts(unitGroup, 1:9);
    h = pie(sizeGroup);
    h = findobj(h, 'Type', 'patch');
    for nh = 1:length(h)
        set(h(nh), 'FaceColor', cmap(nh, :));
    end
    setPrint(8, 6, [PlotDir 'C2SModel/SingleUnitsAnova_' DataSetList(nData).name], 'svg')
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
setPrint(3, 4, [PlotDir 'C2SModel/SingleUnitsAnova_Label'])


close all