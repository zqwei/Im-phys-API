%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cell type categorization
% homogenous cell        : no change of selectivity
% heterogenous cell      : change of selectivity
% no-selective cell      : no selectivity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../Func');
setDir;
load ([TempDatDir 'DataListS2CGP43Model.mat']);


if ~exist([PlotDir 'S2CGP43Model'],'dir')
    mkdir([PlotDir 'S2CGP43Model'])
end

% cmap = cbrewer('qual', 'Set1', 3, 'cubic');
cmap = [0.8000    0.8000    0.8000;
       1.0000    0.6000         0;
       0    0.8000         0];

for nData      = 1:2
    load([TempDatDir DataSetList(nData).name '.mat'])    
    unitGroup = getLogPValueTscoreSpikeTime(nDataSet, DataSetList(nData).params);       
    for nUnit = 1:length(nDataSet)
        nDataSet(nUnit).selectivity = unitGroup(nUnit); %#ok<SAGROW>
    end
    save([TempDatDir DataSetList(nData).name '.mat'], 'nDataSet')
    sizeGroup = histcounts(unitGroup, 0:3);
    disp([sizeGroup(2), sizeGroup(3)])
    figure('Visible', 'off');
    groupNames      = {'Non.', 'Homo.', 'Dynamical'};
    pie(sizeGroup)
    colormap(cmap)
    set(gca, 'TickDir', 'out')
    setPrint(8, 6, [PlotDir 'S2CGP43Model/SingleUnitsTscoreTime_' DataSetList(nData).name])
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cell type categorization
% homogenous cell        : no change of selectivity
% heterogenous cell      : change of selectivity
% no-selective cell      : no selectivity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../Func');
setDir;
load ([TempDatDir 'DataListShuffle.mat']);
nData          = 1;
load([TempDatDir DataSetList(nData).name '.mat'])
unitSpikeGroup = getLogPValueTscoreSpikeTime(nDataSet, DataSetList(nData).params);
spikeDataSet    = nDataSet;
sigma                         = 0.15 / DataSetList(1).params.binsize; % 200 ms
filterLength                  = 10;
filterStep                    = linspace(-filterLength / 2, filterLength / 2, filterLength);
filterInUse                   = exp(-filterStep .^ 2 / (2 * sigma ^ 2));
filterInUse                   = filterInUse / sum (filterInUse);

load ([TempDatDir 'DataListS2CGP43Model.mat']);
if ~exist([PlotDir 'S2CGP43Model'],'dir')
    mkdir([PlotDir 'S2CGP43Model'])
end


% cmap = cbrewer('qual', 'Set1', 3, 'cubic');

for nData      = 1:2
    load([TempDatDir DataSetList(nData).name '.mat'])
    unitGroup = getLogPValueTscoreSpikeTime(nDataSet, DataSetList(nData).params);  
    
    figure; 
    
    groupCounts = zeros(3, 3);
    
    for nGroup  = 1:3
        for mGroup = 1:3
            groupCounts(nGroup, mGroup) = sum(unitGroup == mGroup-1 & unitSpikeGroup == nGroup-1);
        end
    end
    
    groupPerCounts = bsxfun(@rdivide, groupCounts, sum(groupCounts, 2));
    
    figure;
    for nPlot = 1:3
        subplot(1, 3, nPlot)
        bar(1:3, groupPerCounts(nPlot, :), 'edgecolor', 'none');
        ylim([0 1])
        box off
        xlabel('cell type')
        ylabel('frac. cell')
%         colormap(cmap)
        xlim([0.5 3.5])
        set(gca, 'TickDir', 'out')
    end
    setPrint(8*3, 6, [PlotDir 'S2CGP43Model/SingleUnitsTscoreTimeTrans_' DataSetList(nData).name])
        
end

close all