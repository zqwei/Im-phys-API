%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cell type categorization
% homogenous cell        : no change of selectivity
% heterogenous cell      : change of selectivity
% no-selective cell      : no selectivity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../Func');
setDir;
load ([TempDatDir 'DataListC2SRandomDeconvModel.mat']);

% cmap = cbrewer('qual', 'Set1', 3, 'cubic');
cmap = [0.8000    0.8000    0.8000;
       1.0000    0.6000         0;
       0    0.8000         0];   
   
   
   
 
% ffactor     = [0.7; 0.6];
for nData      = 3%1:2
%     load([TempDatDir DataSetList(nData).name '_withOLRemoval.mat'])
    load([TempDatDir DataSetList(nData).name '.mat'])
    fracNeuron = [];
    for ffactor = 0.05:0.05:0.95
        unitGroup = getLogPValueTscoreSpikeTimeAve(nDataSet, DataSetList(nData).params, ffactor);    
        sizeGroup = histcounts(unitGroup, 0:3);
        fracNeuron = [fracNeuron; sizeGroup(2:3)/sum(sizeGroup)];
    end
%     sizeGroup
%     sum(sizeGroup)
%     disp(sizeGroup(2)/sizeGroup(3))
%     figure;
%     groupNames      = {'Non.', 'Homo.', 'Dynamical'};
%     pie(sizeGroup)
%     colormap(cmap)
%     set(gca, 'TickDir', 'out')
%     setPrint(8, 6, [PlotDir 'SingleUnitsTscore/SingleUnitsTscoreTime_' DataSetList(nData).name '_withOLRemoval'])

    figure;
    plot(0.05:0.05:0.95, fracNeuron, '-o')
    xlabel('variability -- Fano factor')
    ylabel('Frac. neuron')
    legend('Mono.', 'Multi.')
    xlim([0 1])
    ylim([0 1])
    set(gca, 'TickDir', 'out')
    setPrint(8, 6, [PlotDir 'SingleUnitsTscore/SingleUnitsTscoreTimeVaribaility_' DataSetList(nData).name '_withOLRemoval'])
end   
   
   

ffactor     = [0.7; 0.6; 0.55];
for nData      = 3 %1:2
%     load([TempDatDir DataSetList(nData).name '_withOLRemoval.mat'])
    load([TempDatDir DataSetList(nData).name '.mat'])
    unitGroup = getLogPValueTscoreSpikeTimeAve(nDataSet, DataSetList(nData).params, ffactor(nData));    
    sizeGroup = histcounts(unitGroup, 0:3);
%     sizeGroup
%     sum(sizeGroup)
%     disp(sizeGroup(2)/sizeGroup(3))
    figure;
    groupNames      = {'Non.', 'Homo.', 'Dynamical'};
    pie(sizeGroup)
    colormap(cmap)
    set(gca, 'TickDir', 'out')
    setPrint(8, 6, [PlotDir 'SingleUnitsTscore/SingleUnitsTscoreTime_' DataSetList(nData).name '_withOLRemoval'])
end

load ([TempDatDir 'DataListShuffle.mat']);
caDataSetList = DataSetList;
load ([TempDatDir 'DataListC2SRandomDeconvModel.mat']);
spikeDataSetList = DataSetList;

ffactor     = [0.7; 0.3];
for nData   = 1:2
    load([TempDatDir caDataSetList(nData + 2).name '_withOLRemoval.mat'])
    caDataSet          = nDataSet;
    unitCaGroup        = getLogPValueTscoreSpikeTime(caDataSet, DataSetList(nData).params); 
    
    load([TempDatDir spikeDataSetList(nData).name '_withOLRemoval.mat'])
    unitGroup = getLogPValueTscoreSpikeTimeAve(nDataSet, DataSetList(nData).params, ffactor(nData)); 
    
    groupCounts = zeros(3, 3);
    for nGroup  = 1:3
        for mGroup = 1:3
            groupCounts(nGroup, mGroup) = sum(unitGroup == mGroup-1 & unitCaGroup == nGroup-1);
        end
    end

    groupPerCounts = bsxfun(@rdivide, groupCounts, sum(groupCounts, 2));

    figure;
    for nPlot = 1:3
        subplot(1, 3, nPlot)
        hold on
        bar(1:3, groupPerCounts(nPlot, :), 'edgecolor', 'none');
        ylim([0 1])
        box off
        xlabel('cell type')
        ylabel('frac. cell')
        colormap(cmap)
        xlim([0.5 3.5])
        set(gca, 'TickDir', 'out')
    end
end

sigma                         = 0.15 / DataSetList(1).params.binsize; % 200 ms
filterLength                  = 10;
filterStep                    = linspace(-filterLength / 2, filterLength / 2, filterLength);
filterInUse                   = exp(-filterStep .^ 2 / (2 * sigma ^ 2));
filterInUse                   = filterInUse / sum (filterInUse);

params = DataSetList(nData).params;
nCell = 227;
figure
subplot(1, 2, 1)
hold on
nUnitData        = caDataSet(nCell).unit_yes_trial;
nUnitData      = getGaussianPSTH (filterInUse, nUnitData, 2);
shadedErrorBar(params.timeSeries, mean(nUnitData, 1),...
    std(nUnitData, 1)/sqrt(size(nUnitData, 1)),...
    {'b-', 'linewid', 1.0}, 0.5);
nUnitData        = caDataSet(nCell).unit_no_trial;
nUnitData      = getGaussianPSTH (filterInUse, nUnitData, 2);
shadedErrorBar(params.timeSeries, mean(nUnitData, 1),...
    std(nUnitData, 1)/sqrt(size(nUnitData, 1)),...
    {'r-', 'linewid', 1.0}, 0.5);

subplot(1, 2, 2)
hold on
nUnitData        = nDataSet(nCell).unit_yes_trial;
nUnitData      = getGaussianPSTH (filterInUse, nUnitData, 2);
shadedErrorBar(params.timeSeries, mean(nUnitData, 1),...
    std(nUnitData, 1)/sqrt(size(nUnitData, 1)),...
    {'b-', 'linewid', 1.0}, 0.5);
nUnitData        = nDataSet(nCell).unit_no_trial;
nUnitData      = getGaussianPSTH (filterInUse, nUnitData, 2);
shadedErrorBar(params.timeSeries, mean(nUnitData, 1),...
    std(nUnitData, 1)/sqrt(size(nUnitData, 1)),...
    {'r-', 'linewid', 1.0}, 0.5);
