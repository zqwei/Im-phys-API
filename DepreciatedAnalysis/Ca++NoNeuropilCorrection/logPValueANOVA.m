%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cell type categorization
% AP cell                : sensory
% LR cell                : decision
% CE cell                : reward
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../Func');
setDir;
load ([TempDatDir 'DataListNoNeuropil.mat']);
cmap = cbrewer('qual', 'Set1', 10, 'cubic');

for nData      = 1:length(DataSetList)
    load([TempDatDir DataSetList(nData).name '.mat'])
    logPValueEpoch= getLogPValueSpikeEpoch(nDataSet, DataSetList(nData).params);
    unitGroup     = getAnovaLogPSpikeEpoch(logPValueEpoch);
    figure
    groupNames    = {'Pole', 'Lick', 'Reward', 'PL', 'PR', 'LR', 'PLR', 'Non.'};
    sizeGroup     = histcounts(unitGroup, 1:9);
    h = pie(sizeGroup);
    h = findobj(h, 'Type', 'patch');
    for nh = 1:length(h)
        set(h(nh), 'FaceColor', cmap(nh, :));
    end
    setPrint(8, 6, [PlotDir 'NoNeuropil/SingleUnitsAnova_' DataSetList(nData).name], 'svg')

end

close all