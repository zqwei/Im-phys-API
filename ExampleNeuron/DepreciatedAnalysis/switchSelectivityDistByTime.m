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


nData   = 1;
load('S2CC2S.mat', 'nDataSet')

for nUnit        = 1:length(nDataSet)
    nDataSet(nUnit).unit_yes_trial      = nDataSet(nUnit).mcmc_yes_trial(:, 56:132);
    nDataSet(nUnit).unit_no_trial       = nDataSet(nUnit).mcmc_no_trial(:, 56:132);
end
unitGroup = getLogPValueTscoreSpikeTime(nDataSet, DataSetList(nData).params); 
sizeGroup = histcounts(unitGroup, 0:3);

figure('Visible', 'off');
groupNames      = {'Non.', 'Homo.', 'Dynamical'};
pie(sizeGroup)
% colormap(cmap)
set(gca, 'TickDir', 'out')
setPrint(8, 6, 'SingleUnitsTscoreTime_MCMC')


addpath('../Func');
setDir;
load([TempDatDir 'Shuffle_Spikes.mat'])
unitSpikeGroup = getLogPValueTscoreSpikeTime(nDataSet, DataSetList(nData).params);
spikeDataSet    = nDataSet;
sigma                         = 0.15 / DataSetList(1).params.binsize; % 200 ms
filterLength                  = 10;
filterStep                    = linspace(-filterLength / 2, filterLength / 2, filterLength);
filterInUse                   = exp(-filterStep .^ 2 / (2 * sigma ^ 2));
filterInUse                   = filterInUse / sum (filterInUse);
load ([TempDatDir 'DataListS2CModel.mat']);

nData   = 1;

load('S2CC2S.mat', 'nDataSet')
    
for nUnit        = 1:length(nDataSet)
    nDataSet(nUnit).unit_yes_trial      = nDataSet(nUnit).mcmc_yes_trial(:, 56:132);
    nDataSet(nUnit).unit_no_trial       = nDataSet(nUnit).mcmc_no_trial(:, 56:132);
end
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
setPrint(8*3, 6, ['SingleUnitsTscoreTimeTrans_S2CC2S_' DataSetList(nData).name])
