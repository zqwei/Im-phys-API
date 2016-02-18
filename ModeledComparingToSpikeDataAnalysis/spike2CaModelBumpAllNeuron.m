%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% low firing activity to calcium analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% function spike2CaModelBumpAllNeuron


Fm              = 21.3240;
K               = 13.9248;
n               = 1.7531;
tau_rise        = 0.0728;
tau_decay       = 1.4551;
intNoise        = 2.5;
extNoise        = 0;
nFactor         = 0.5:0.1:1.5;

addpath('../Func');
setDir;
load ([TempDatDir 'DataListModeled.mat']);
nData = 1;
load([TempDatDir DataSetList(nData).name '.mat'])
params                 = DataSetList(nData).params;
params.Fm              = Fm;
params.K               = K;
params.n               = n;

SSB_SST_mat     = nan(length(nDataSet), 3);

for nNeuron     = 1:length(nDataSet)
    spikeDataSet           = nDataSet(nNeuron);
    
    nYesTrial              = nan(length(nFactor), size(spikeDataSet.unit_yes_trial, 2));
    nNoTrial               = nan(length(nFactor), size(spikeDataSet.unit_yes_trial, 2));

    for nTau_rise              = 1:length(nFactor)
        params.tau_r           = tau_rise * nFactor(nTau_rise);
        params.tau_d           = tau_decay;
        params.intNoise        = intNoise;
        params.extNoise        = extNoise;
        caDataSet               = getFakeCaImagingData(spikeDataSet, params);
        nYesTrial(nTau_rise,:) = mean(caDataSet.unit_yes_trial);
        nNoTrial(nTau_rise,:)  = mean(caDataSet.unit_no_trial);
    end
    
    SSB_SST_mat(nNeuron, 1)    = (mean(var(nYesTrial))/var(nYesTrial(:)) + mean(var(nNoTrial))/var(nNoTrial(:)))/2;

    nYesTrial              = nan(length(nFactor), size(spikeDataSet.unit_yes_trial, 2));
    nNoTrial               = nan(length(nFactor), size(spikeDataSet.unit_yes_trial, 2));
    for nTau_decay              = 1:length(nFactor)
        params.tau_r           = tau_rise;
        params.tau_d           = tau_decay * nFactor(nTau_decay);
        params.intNoise        = intNoise;
        params.extNoise        = extNoise;
        caDataSet               = getFakeCaImagingData(spikeDataSet, params);
        nYesTrial(nTau_decay,:) = mean(caDataSet.unit_yes_trial);
        nNoTrial(nTau_decay,:)  = mean(caDataSet.unit_no_trial);
    end
    
    SSB_SST_mat(nNeuron, 2)    = (mean(var(nYesTrial))/var(nYesTrial(:)) + mean(var(nNoTrial))/var(nNoTrial(:)))/2;

    nYesTrial              = nan(length(nFactor), size(spikeDataSet.unit_yes_trial, 2));
    nNoTrial               = nan(length(nFactor), size(spikeDataSet.unit_yes_trial, 2));
    for nintNoise              = 1:length(nFactor)
        params.tau_r           = tau_rise;
        params.tau_d           = tau_decay;
        params.intNoise        = intNoise * nFactor(nintNoise);
        params.extNoise        = extNoise;
        caDataSet               = getFakeCaImagingData(spikeDataSet, params);
        nYesTrial(nintNoise,:) = mean(caDataSet.unit_yes_trial);
        nNoTrial(nintNoise,:)  = mean(caDataSet.unit_no_trial);
    end
    
    SSB_SST_mat(nNeuron, 3)    = (mean(var(nYesTrial))/var(nYesTrial(:)) + mean(var(nNoTrial))/var(nNoTrial(:)))/2;
end

figure;
subplot(1, 2, 1)
plot(SSB_SST_mat(:, 1), SSB_SST_mat(:, 2), 'ok', 'linewid', 1)
h = refline([1 0]);
h.Color = 'r';
h.LineStyle = '--';
h.LineWidth = 2;
xlabel('% SSB variance (\tau_{rise})')
ylabel('% SSB variance (\tau_{decay})')
set(gca, 'TickDir', 'out')
box off


subplot(1, 2, 2)
plot(SSB_SST_mat(:, 1), SSB_SST_mat(:, 3), 'ok', 'linewid', 1)
h = refline([1 0]);
h.Color = 'r';
h.LineStyle = '--';
h.LineWidth = 2;
xlabel('% SSB variance (\tau_{rise})')
ylabel('% SSB variance (\tau_{decay})')
set(gca, 'TickDir', 'out')
box off

setPrint(8*2, 6, 'Spike2CaAnalysisPerSSB', 'pdf')