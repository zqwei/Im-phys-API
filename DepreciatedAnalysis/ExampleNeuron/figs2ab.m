addpath('../Func');
setDir;

load ([TempDatDir 'DataListShuffle.mat']);
nData = 1; % plot raster and psth
load([TempDatDir DataSetList(nData).name '.mat'])
spikeDataSet = nDataSet;   
params               = DataSetList(nData).params;

sigma                         = 0.15 / params.binsize; % 200 ms
filterLength                  = 11;
filterStep                    = linspace(-filterLength / 2, filterLength / 2, filterLength);
filterInUse                   = exp(-filterStep .^ 2 / (2 * sigma ^ 2));
filterInUse                   = filterInUse / sum (filterInUse); 

numTrials            = 50;
spkDataSet           = nDataSet;
spkFiringRates       = generateDPCAData(spkDataSet, numTrials);
spkFiringRatesMean   = squeeze(mean(spkFiringRates, 4));
spkFiringRatesMean   = [squeeze(spkFiringRatesMean(:, 1, :)), squeeze(spkFiringRatesMean(:, 2, :))];



load([TempDatDir 'FineTunedModeled_6s_AAV.mat']);
s2cDataSet           = nDataSet;
s2cFiringRates       = generateDPCAData(s2cDataSet, numTrials);
s2cFiringRatesMean   = squeeze(mean(s2cFiringRates, 4));
s2cFiringRatesMean   = [squeeze(s2cFiringRatesMean(:, 1, :)), squeeze(s2cFiringRatesMean(:, 2, :))];

rhoSpk               = corr(spkFiringRatesMean', spkFiringRatesMean', 'type', 'Spearman');
rhoS2c               = corr(s2cFiringRatesMean', s2cFiringRatesMean', 'type', 'Spearman');

neuronIndex          = [375 46];
for nCell            = 1:length(neuronIndex)
    nCellid          = neuronIndex(nCell);
    subplot(length(neuronIndex), 2, 2*(nCell-1)+1)
    hold on;
    nUnitData        = spikeDataSet(nCellid).unit_yes_trial;
    nUnitData        = getGaussianPSTH (filterInUse, nUnitData, 2);
    shadedErrorBar(params.timeSeries, mean(nUnitData, 1),...
        std(nUnitData, 1)/sqrt(size(nUnitData, 1)),...
        {'k-', 'linewid', 1.0, 'color', 'b'}, 0.5);
    nUnitData        = spikeDataSet(nCellid).unit_no_trial;
    nUnitData        = getGaussianPSTH (filterInUse, nUnitData, 2);
    shadedErrorBar(params.timeSeries, mean(nUnitData, 1),...
        std(nUnitData, 1)/sqrt(size(nUnitData, 1)),...
        {'k-', 'linewid', 1.0, 'color', 'r'}, 0.5);
    gridxy ([params.polein, params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 1.0)
    hold off;
    ylabel('Spikes /s');
    xlabel('Time (s)');
    xlim([params.timeSeries(1) params.timeSeries(end)]);
    set(gca, 'TickDir', 'out')
    subplot(length(neuronIndex), 2, 2*(nCell-1)+2)
    hold on;
    nUnitData        = s2cDataSet(nCellid).unit_yes_trial;
    nUnitData        = getGaussianPSTH (filterInUse, nUnitData, 2);
    shadedErrorBar(params.timeSeries, mean(nUnitData, 1),...
        std(nUnitData, 1)/sqrt(size(nUnitData, 1)),...
        {'k-', 'linewid', 1.0, 'color', 'b'}, 0.5);
    nUnitData        = s2cDataSet(nCellid).unit_no_trial;
    nUnitData        = getGaussianPSTH (filterInUse, nUnitData, 2);
    shadedErrorBar(params.timeSeries, mean(nUnitData, 1),...
        std(nUnitData, 1)/sqrt(size(nUnitData, 1)),...
        {'k-', 'linewid', 1.0, 'color', 'r'}, 0.5);
    gridxy ([params.polein, params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 1.0)
    hold off;
    ylabel('df/f');
    xlabel('Time (s)');
    xlim([params.timeSeries(1) params.timeSeries(end)]);
    set(gca, 'TickDir', 'out')
end