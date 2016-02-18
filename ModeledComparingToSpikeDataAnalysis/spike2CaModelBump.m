%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% low firing activity to calcium analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function spike2CaModelBump(nNeuron)


Fm              = 21.3240;
K               = 13.9248;
n               = 1.7531;
tau_rise        = 0.0728;
tau_decay       = 1.4551;
intNoise        = 2.5;
extNoise        = 0;
nFactor         = 0.5:0.1:1.5;
colorSet        = cool(length(nFactor));

addpath('../Func');
setDir;
load ([TempDatDir 'DataListModeled.mat']);
nData = 1;
load([TempDatDir DataSetList(nData).name '.mat'])

%%%%%%% plot discrete-time-binned psth
% mkdir('SpikeActPlots')
% 
% for nNeuron = 1:length(nDataSet)
%     h = figure;
%     plotMeanActivityTrace (nDataSet, nNeuron, DataSetList(nData).params, 'Firing Rate (Hz)', 'Time (s)')
%     setPrint(8, 6, ['SpikeActPlots/SpikeActPlots_' num2str(nNeuron, '%04d')], 'tif')
%     close(h)
% end

%%%% pick a neuronal index in the population

% nNeuron                = 17;
spikeDataSet           = nDataSet(nNeuron); %#ok<NODEF>
params                 = DataSetList(nData).params;
params.Fm              = Fm;
params.K               = K;
params.n               = n;


% plot spiking raster
figure;
spkTimes{1}    = spikeDataSet.unit_yes_trial_spk_time;
spkTimes{2}    = spikeDataSet.unit_no_trial_spk_time;
color_index    = [0.7 0.0 0.0; 0.0 0.0 0.7]; % CL; CR

subplot(3, 3, [1 4])
spkLoc         = 0;
for nPlot            = 1:2
    hold on;
    for nTrial       = 1:length(spkTimes{nPlot})
        spikeTimes   = spkTimes{nPlot}{nTrial};
        spikeTrial   = ones(length(spikeTimes), 1) * (nTrial + spkLoc);
        plot(spikeTimes, spikeTrial, '.', 'color', color_index(nPlot, :));
    end
    spkLoc     = spkLoc + length(spkTimes{nPlot}) + 3;
end

gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 1.0);
hold off;
ylim([1 spkLoc])
xlim([DataSetList(nData).params.timeSeries(1) DataSetList(nData).params.timeSeries(end)]);
axis off

% plot spiking psth
subplot(3, 3, 7)
hold on;
nUnitData        = spikeDataSet.unit_yes_trial;
shadedErrorBar(DataSetList(nData).params.timeSeries, mean(nUnitData, 1),...
    std(nUnitData, 1)/sqrt(size(nUnitData, 1)),...
    {'k-', 'linewid', 1.0, 'color', color_index(1, :)}, 0.5);
nUnitData        = spikeDataSet.unit_no_trial;
shadedErrorBar(DataSetList(nData).params.timeSeries, mean(nUnitData, 1),...
    std(nUnitData, 1)/sqrt(size(nUnitData, 1)),...
    {'k-', 'linewid', 1.0, 'color', color_index(2, :)}, 0.5);
gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 1.0)
hold off;
ylabel('Firing rate (Hz)');
xlabel('Time (s)');
xlim([DataSetList(nData).params.timeSeries(1) DataSetList(nData).params.timeSeries(end)]);




% change of tau_rise



for nTau_rise              = 1:length(nFactor)
    params.tau_r           = tau_rise * nFactor(nTau_rise);
    params.tau_d           = tau_decay;
    params.intNoise        = intNoise;
    params.extNoise        = extNoise;
    nDataSet               = getFakeCaImagingData(spikeDataSet, params);
    
    subplot(3, 3, 2)
    hold on;
    nUnitData        = nDataSet.unit_yes_trial;
    plot(DataSetList(nData).params.timeSeries, mean(nUnitData, 1),'-', 'linewid', 1.0, 'color', colorSet(nTau_rise, :));
    gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 1.0)
    hold off;
    ylabel('DF/F (\tau_{rise})');
    xlabel('Time (s)');
    xlim([DataSetList(nData).params.timeSeries(1) DataSetList(nData).params.timeSeries(end)]);
    
    subplot(3, 3, 3);
    hold on;
    nUnitData        = nDataSet.unit_no_trial;
    plot(DataSetList(nData).params.timeSeries, mean(nUnitData, 1),'-', 'linewid', 1.0, 'color', colorSet(nTau_rise, :));
    gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 1.0)
    hold off;
    ylabel('DF/F (\tau_{rise})');
    xlabel('Time (s)');
    xlim([DataSetList(nData).params.timeSeries(1) DataSetList(nData).params.timeSeries(end)]);
end


for nTau_decay              = 1:length(nFactor)
    params.tau_r           = tau_rise;
    params.tau_d           = tau_decay * nFactor(nTau_decay);
    params.intNoise        = intNoise;
    params.extNoise        = extNoise;
    nDataSet               = getFakeCaImagingData(spikeDataSet, params);
    
    subplot(3, 3, 5)
    hold on;
    nUnitData        = nDataSet.unit_yes_trial;
    plot(DataSetList(nData).params.timeSeries, mean(nUnitData, 1),'-', 'linewid', 1.0, 'color', colorSet(nTau_decay, :));
    gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 1.0)
    hold off;
    ylabel('DF/F (\tau_{decay})');
    xlabel('Time (s)');
    xlim([DataSetList(nData).params.timeSeries(1) DataSetList(nData).params.timeSeries(end)]);
    
    subplot(3, 3, 6)
    hold on;
    nUnitData        = nDataSet.unit_no_trial;
    plot(DataSetList(nData).params.timeSeries, mean(nUnitData, 1),'-', 'linewid', 1.0, 'color', colorSet(nTau_decay, :));
    gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 1.0)
    hold off;
    ylabel('DF/F (\tau_{decay})');
    xlabel('Time (s)');
    xlim([DataSetList(nData).params.timeSeries(1) DataSetList(nData).params.timeSeries(end)]);
end


for nintNoise              = 1:length(nFactor)
    params.tau_r           = tau_rise;
    params.tau_d           = tau_decay;
    params.intNoise        = intNoise * nFactor(nintNoise);
    params.extNoise        = extNoise;
    nDataSet               = getFakeCaImagingData(spikeDataSet, params);
    
    subplot(3, 3, 8)
    hold on;
    nUnitData        = nDataSet.unit_yes_trial;
    plot(DataSetList(nData).params.timeSeries, mean(nUnitData, 1),'-', 'linewid', 1.0, 'color', colorSet(nintNoise, :));
    gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 1.0)
    hold off;
    ylabel('DF/F (\sigma_{int})');
    xlabel('Time (s)');
    xlim([DataSetList(nData).params.timeSeries(1) DataSetList(nData).params.timeSeries(end)]);
    
    subplot(3, 3, 9)
    hold on;
    nUnitData        = nDataSet.unit_no_trial;
    plot(DataSetList(nData).params.timeSeries, mean(nUnitData, 1),'-', 'linewid', 1.0, 'color', colorSet(nintNoise, :));
    gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 1.0)
    hold off;
    ylabel('DF/F (\sigma_{int})');
    xlabel('Time (s)');
    xlim([DataSetList(nData).params.timeSeries(1) DataSetList(nData).params.timeSeries(end)]);
end

setPrint(8*3, 6*3, ['Spike2CaAnalysis_' num2str(nNeuron, '%04d')], 'pdf')