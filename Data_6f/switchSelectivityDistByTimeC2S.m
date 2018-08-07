%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cell type categorization
% homogenous cell        : no change of selectivity
% heterogenous cell      : change of selectivity
% no-selective cell      : no selectivity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../Func');
setDir;
load([TempDatDir 'DataListShuffle.mat']);
load([TempDatDir 'Shuffle_Spikes_Nuo_Short_Delay.mat'])
nData          = 1;
depth          = [nDataSet.depth_in_um];
spikeDataSet   = nDataSet(depth<471);
unitSpikeGroup = getLogPValueTscoreSpikeTime(spikeDataSet, DataSetList(nData).params);  
sigma                         = 0.15 / DataSetList(1).params.binsize; % 200 ms
filterLength                  = 10;
filterStep                    = linspace(-filterLength / 2, filterLength / 2, filterLength);
filterInUse                   = exp(-filterStep .^ 2 / (2 * sigma ^ 2));
filterInUse                   = filterInUse / sum (filterInUse);

fName = 'Ca_Fast_SShort_Delay';

load(['S2CC2S_Deconv_' fName '.mat'])

spkDataSet = nDataSet;

% fracNeuron = [];
% for ffactor = 0.5:0.05:3.95
%     unitGroup = getLogPValueTscoreSpikeTimeAve(spkDataSet, DataSetList(nData).params, ffactor);    
%     sizeGroup = histcounts(unitGroup, 0:3);
%     fracNeuron = [fracNeuron; sizeGroup(2:3)/sum(sizeGroup)];
% end
% 
% figure;
% plot(0.5:0.05:3.95, fracNeuron, '-o')
% xlabel('variability -- Fano factor')
% ylabel('Frac. neuron')
% legend('Mono.', 'Multi.')
% xlim([0.49 4])
% ylim([0 1])
% set(gca, 'TickDir', 'out')

ffactor   = 2.5;
unitGroup = getLogPValueTscoreSpikeTimeAve(spkDataSet, DataSetList(nData).params, ffactor);  

figure; 
false_pos = unitGroup ~= unitSpikeGroup;
unitGroup = unitGroup*2;
unitGroup(false_pos) = unitGroup(false_pos)+1;
sizeGroup = histcounts(unitGroup, 0:6);
pie(sizeGroup)
disp(sizeGroup/sum(sizeGroup))
setPrint(8, 6, [PlotDir 'SingleUnitsTscore/SingleUnitsTscoreTime_S2CC2S_' fName])
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
    xlim([0.5 3.5])
    set(gca, 'TickDir', 'out')
end
setPrint(8*3, 6, [PlotDir 'SingleUnitsTscore/SingleUnitsTscoreTimeTrans_S2CC2S_' fName])