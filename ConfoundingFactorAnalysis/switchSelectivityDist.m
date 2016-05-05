%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cell type categorization
% homogenous cell        : no change of selectivity
% heterogenous cell      : change of selectivity
% no-selective cell      : no selectivity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../Func');
setDir;
load ([TempDatDir 'DataListShuffleConfounding.mat']);


if ~exist([PlotDir 'ConfoundingFactorTscore'],'dir')
    mkdir([PlotDir 'ConfoundingFactorTscore'])
end

cmap = cbrewer('qual', 'Set1', 3, 'cubic');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All area vs area of spiking recording
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nData          = 5;
load([TempDatDir 'LogPValueTscore_' DataSetList(nData).name '.mat'], 'logPValueEpoch')
unitGroup      = plotTtestLogPSpikeEpoch (logPValueEpoch);
figure;
sizeGroup      = histcounts(unitGroup, 0:3);
pie(sizeGroup)
colormap(cmap)
disp(sum(sizeGroup))
setPrint(8, 6, [PlotDir 'ConfoundingFactorTscore/SingleUnitsTscoreAPML_Pop_' DataSetList(nData).name])

figure;
APLoc          = [DataSetList(nData).cellinfo.AP_axis];
MLLoc          = [DataSetList(nData).cellinfo.ML_axis];
sizeGroup      = histcounts(unitGroup(APLoc>2400 & APLoc<2600 & MLLoc>1100 & MLLoc<1900), 0:3);
groupNames     = {'Non.', 'Homo.', 'Dynamicial'};
pie(sizeGroup)
colormap(cmap)
disp(sum(sizeGroup))
setPrint(8, 6, [PlotDir 'ConfoundingFactorTscore/SingleUnitsTscoreAPML_Sub_' DataSetList(nData).name])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison across animals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load ([TempDatDir 'DataListShuffleConfounding.mat']);
nData      = 1;    
[~, ~, anmIndex] = unique(cell2mat({DataSetList(nData).cellinfo.anmName}'), 'rows');
numGroup    = 3;
groupCounts = zeros(anmIndex(end), numGroup);
load([TempDatDir 'LogPValueTscore_' DataSetList(nData).name '.mat'], 'logPValueEpoch')
unitGroup   = plotTtestLogPSpikeEpoch (logPValueEpoch);
for nAnm    = 1:anmIndex(end)
    if sum(anmIndex == nAnm) > 50
        for nGroup = 1:numGroup
            groupCounts(nAnm, nGroup) = sum(unitGroup == nGroup-1 & anmIndex == nAnm);
        end
    end
end
zeroGroups  = sum(groupCounts, 2) == 0;
groupCounts = groupCounts(~zeroGroups, :);
groupPerCounts = bsxfun(@rdivide, groupCounts, sum(groupCounts, 2));
figure
subplot(2, 1, 1)
barh(groupPerCounts, 'stack', 'edgecolor', 'k');
colormap(cmap)
%     caxis([1 8])
xlim([0 1])
box off
xlabel('% cell type')
ylabel('Animal index')
set(gca, 'yTickLabel', {})
subplot(2, 1, 2)
barh(sum(groupCounts,2), 'k')
xlabel('# cells')
ylabel('Animal index')
set(gca, 'yTickLabel', {})
setPrint(8, 6*2, [PlotDir 'ConfoundingFactorTscore/SingleUnitsTscoreANM_' DataSetList(nData).name])


nData       = 6;    
[~, ~, anmIndex] = unique(cell2mat({DataSetList(nData).cellinfo.anmName}'), 'rows');
load ([TempDatDir 'DataListShuffle.mat']);
nData       = 4;
load([TempDatDir DataSetList(nData).name '_withOLRemoval.mat'])
logPValueEpoch= getLogPValueTscoreSpikeEpoch(nDataSet, DataSetList(nData).params);
anmIndex    = anmIndex(~neuronRemoveList);
numGroup    = 3;
groupCounts = zeros(anmIndex(end), numGroup);
unitGroup   = plotTtestLogPSpikeEpoch (logPValueEpoch);
for nAnm    = 1:anmIndex(end)
    if sum(anmIndex == nAnm) > 50
        for nGroup = 1:numGroup
            groupCounts(nAnm, nGroup) = sum(unitGroup == nGroup-1 & anmIndex == nAnm);
        end
    end
end

zeroGroups  = sum(groupCounts, 2) == 0;
groupCounts = groupCounts(~zeroGroups, :);
groupPerCounts = bsxfun(@rdivide, groupCounts, sum(groupCounts, 2));
figure
subplot(2, 1, 1)
barh(groupPerCounts, 'stack', 'edgecolor', 'k');
colormap(cmap)
%     caxis([1 8])
xlim([0 1])
box off
xlabel('% cell type')
ylabel('Animal index')
set(gca, 'yTickLabel', {})

subplot(2, 1, 2)
barh(sum(groupCounts,2), 'k')
xlabel('# cells')
ylabel('Animal index')
set(gca, 'yTickLabel', {})

setPrint(8, 6*2, [PlotDir 'ConfoundingFactorTscore/SingleUnitsTscoreANM_' DataSetList(nData).name])

close all
