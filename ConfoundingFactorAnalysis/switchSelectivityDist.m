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


DataSetToAnalysis = [1 5 6];

% for nData      = DataSetToAnalysis
%     load([TempDatDir DataSetList(nData).name '.mat'])
%     logPValueEpoch= getLogPValueTscoreSpikeEpoch(nDataSet, DataSetList(nData).params);
%     save([TempDatDir 'LogPValueTscore_' DataSetList(nData).name '.mat'], 'logPValueEpoch')
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All area vs area of spiking recording
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nData          = 5;
load([TempDatDir 'LogPValueTscore_' DataSetList(nData).name '.mat'], 'logPValueEpoch')
unitGroup      = plotTtestLogPSpikeEpoch (logPValueEpoch);
figure;
subplot(1, 2, 1)
sizeGroup      = histcounts(unitGroup, 0:3);
groupNames     = {'Non.', 'Homo.', 'Dynamicial'};
pie(sizeGroup, groupNames)
title('Population')

APLoc          = [DataSetList(nData).cellinfo.AP_axis];
MLLoc          = [DataSetList(nData).cellinfo.ML_axis];
subplot(1, 2, 2)
sizeGroup      = histcounts(unitGroup(APLoc>2400 & APLoc<2600 & MLLoc>1100 & MLLoc<1900), 0:3);
groupNames     = {'Non.', 'Homo.', 'Dynamicial'};
pie(sizeGroup, groupNames)
title('Subpopulation')

setPrint(8*2, 6, [PlotDir 'ConfoundingFactorTscore/SingleUnitsTscoreAPML_' DataSetList(nData).name], 'pdf')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison across animals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for nData      = DataSetToAnalysis    
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
    subplot(1, 2, 1)
    barh(groupPerCounts, 'stack', 'edgecolor', 'none');
%     caxis([1 8])
    xlim([0 1])
    box off
    xlabel('% cell type')
    ylabel('Animal index')
    set(gca, 'yTickLabel', {})

    subplot(1, 2, 2)
    barh(sum(groupCounts,2))
    xlabel('# cells')
    ylabel('Animal index')
    set(gca, 'yTickLabel', {})
    
    setPrint(8*2, 6, [PlotDir 'ConfoundingFactorTscore/SingleUnitsTscoreANM_' DataSetList(nData).name], 'pdf')
    
end

close all
