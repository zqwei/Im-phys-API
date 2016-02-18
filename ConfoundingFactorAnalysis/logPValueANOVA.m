%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cell type categorization
% AP cell                : sensory
% LR cell                : decision
% CE cell                : reward
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../Func');
setDir;
load ([TempDatDir 'DataListShuffleConfounding.mat']);


if ~exist([PlotDir 'ConfoundingFactorAnova'],'dir')
    mkdir([PlotDir 'ConfoundingFactorAnova'])
end

DataSetToAnalysis = [1 5 6];

% for nData      = DataSetToAnalysis
%     load([TempDatDir DataSetList(nData).name '.mat'])
%     logPValueEpoch= getLogPValueSpikeEpoch(nDataSet, DataSetList(nData).params);
%     save([TempDatDir 'LogPValue_' DataSetList(nData).name '.mat'], 'logPValueEpoch')
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All area vs area of spiking recording
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


nData          = 5;
load([TempDatDir 'LogPValue_' DataSetList(nData).name '.mat'], 'logPValueEpoch')
unitGroup      = groupAnovaLogPSpikeEpoch (logPValueEpoch);
figure;
subplot(1, 2, 1)
sizeGroup      = histcounts(unitGroup, 1:9);
groupNames     = {'Pole', 'Lick', 'Reward', 'PL', 'PR', 'LR', 'PLR', 'Non-Selective'};
pie(sizeGroup, groupNames)
title('Population')

APLoc          = [DataSetList(nData).cellinfo.AP_axis];
MLLoc          = [DataSetList(nData).cellinfo.ML_axis];
subplot(1, 2, 2)
sizeGroup      = histcounts(unitGroup(APLoc>2400 & APLoc<2600 & MLLoc>1100 & MLLoc<1900), 1:9);
groupNames     = {'Pole', 'Lick', 'Reward', 'PL', 'PR', 'LR', 'PLR', 'Non-Selective'};
pie(sizeGroup, groupNames)
title('Subpopulation')

setPrint(8*2, 6, [PlotDir 'ConfoundingFactorAnova/SingleUnitsAnovaAPML_' DataSetList(nData).name], 'pdf')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison across animals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for nData      = DataSetToAnalysis    
    [~, ~, anmIndex] = unique(cell2mat({DataSetList(nData).cellinfo.anmName}'), 'rows');
    numGroup    = 7;
    groupCounts = zeros(anmIndex(end), numGroup);
    load([TempDatDir 'LogPValue_' DataSetList(nData).name '.mat'], 'logPValueEpoch')
    unitGroup   = groupAnovaLogPSpikeEpoch (logPValueEpoch);
    for nAnm    = 1:anmIndex(end)
        if sum(anmIndex == nAnm) > 50
            for nGroup = 1:numGroup
                groupCounts(nAnm, nGroup) = sum(unitGroup == nGroup & anmIndex == nAnm);
            end
        end
    end
    
    zeroGroups  = sum(groupCounts, 2) == 0;
    groupCounts = groupCounts(~zeroGroups, :);
    groupPerCounts = bsxfun(@rdivide, groupCounts, sum(groupCounts, 2));
    figure
    subplot(1, 2, 1)
    barh(groupPerCounts, 'stack', 'edgecolor', 'none');
    caxis([1 8])
    xlim([0 1])
    box off
    xlabel('% cell type')
    ylabel('Animal index')
    set(gca, 'yTickLabel', {})

    subplot(1, 2, 2)
    bar(sum(groupCounts,2))
    ylabel('# cells')
    xlabel('Animal index')
    set(gca, 'xTickLabel', {})
    
    setPrint(8*2, 6, [PlotDir 'ConfoundingFactorAnova/SingleUnitsAnovaANM_' DataSetList(nData).name], 'pdf')
    
end

close all