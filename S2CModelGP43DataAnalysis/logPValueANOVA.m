%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cell type categorization
% AP cell                : sensory
% LR cell                : decision
% CE cell                : reward
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../Func');
setDir;
load ([TempDatDir 'DataListS2CGP43Model.mat']);


cmap = cbrewer('qual', 'Set1', 10, 'cubic');

for nData      = [1 2]
    load([TempDatDir DataSetList(nData).name '.mat'])
    logPValueEpoch= getLogPValueSpikeEpoch(nDataSet, DataSetList(nData).params);
    unitGroup = getAnovaLogPSpikeEpoch(logPValueEpoch);
    figure
    groupNames      = {'Stimulus', 'Movement', 'Outcome', 'Mixed', 'Non.'};
    sizeGroup      = histcounts(unitGroup, 1:5);
    h = pie(sizeGroup);
    h = findobj(h, 'Type', 'patch');
    for nh = 1:length(h)
        set(h(nh), 'FaceColor', cmap(nh, :));
    end
    setPrint(8, 6, [PlotDir 'SingleUnitsAnova/SingleUnitsAnova_' DataSetList(nData).name])

end

close all