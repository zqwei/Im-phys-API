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

cmap = cbrewer('qual', 'Set1', 10, 'cubic');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All area vs area of spiking recording
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


nData          = 5;
load([TempDatDir 'LogPValue_' DataSetList(nData).name '.mat'], 'logPValueEpoch')
unitGroup      = getAnovaLogPSpikeEpoch (logPValueEpoch);
figure;
sizeGroup      = histcounts(unitGroup, 1:9);
h = pie(sizeGroup);
h = findobj(h, 'Type', 'patch');
for nh = 1:length(h)
    set(h(nh), 'FaceColor', cmap(nh, :));
end
setPrint(8, 6, [PlotDir 'ConfoundingFactorAnova/SingleUnitsAnovaAPML_Pop_' DataSetList(nData).name], 'svg')

APLoc          = [DataSetList(nData).cellinfo.AP_axis];
MLLoc          = [DataSetList(nData).cellinfo.ML_axis];
figure;
sizeGroup      = histcounts(unitGroup(APLoc>2400 & APLoc<2600 & MLLoc>1100 & MLLoc<1900), 1:9);
h = pie(sizeGroup);
h = findobj(h, 'Type', 'patch');
for nh = 1:length(h)
    set(h(nh), 'FaceColor', cmap(nh, :));
end
setPrint(8, 6, [PlotDir 'ConfoundingFactorAnova/SingleUnitsAnovaAPML_Sub_' DataSetList(nData).name], 'svg')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison across animals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DataSetToAnalysis = [1 6];
for nData      = DataSetToAnalysis    
    [~, ~, anmIndex] = unique(cell2mat({DataSetList(nData).cellinfo.anmName}'), 'rows');
    numGroup    = 7;
    groupCounts = zeros(anmIndex(end), numGroup);
    load([TempDatDir 'LogPValue_' DataSetList(nData).name '.mat'], 'logPValueEpoch')
    unitGroup   = getAnovaLogPSpikeEpoch (logPValueEpoch);
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
    colormap(cmap(1:8, :))

    subplot(1, 2, 2)
    barh(sum(groupCounts,2), 'k')
    xlabel('# cells')
    ylabel('Animal index')
    set(gca, 'yTickLabel', {})
    
    setPrint(8*2, 6, [PlotDir 'ConfoundingFactorAnova/SingleUnitsAnovaANM_' DataSetList(nData).name], 'svg')
    
end

close all