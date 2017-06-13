%
% Compute distribution of selectivity for each single neuron
% 
% Here we drop the test of parameter of tau_r
% 
% 
% -------------------------------------------------------------------------
% version 1.0
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Common parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('../Func');
setDir;
% % load ([TempDatDir 'DataListShuffle.mat']);
% % numFold = 30;
% % 
% % % short Ca slow
% % nData  = 3;
% % params = DataSetList(nData).params;
% % std_r                  = 0.0375; %0.0375;
% % median_r               = 0.0927;
% % std_d                  = 0.5374; %0.5374;
% % median_d               = 1.2294;
% % load([TempDatDir DataSetList(nData).name '_withOLRemoval.mat'])
% % per_list               = 0.02:0.01:0.98;
% % tau_d_list             = icdf('Normal', per_list, 0, 1) * std_d + median_d;
% % noise_factor_list      = sqrt(icdf('Exponential', per_list, 0.3527));
% % numNeuron              = length(nDataSet);
% % cellTypeMat            = nan(numNeuron, numFold, length(noise_factor_list), length(tau_d_list));
% % for nTau              = 1:length(tau_d_list)
% %     spikeDataSet      = getSyntheticSpikeDeconvDataSimpleVersion(nDataSet, median_r, tau_d_list(nTau), params);  
% %     for nFF           = 1:length(noise_factor_list) 
% %         for nFold     = 1:numFold
% %             unitGroup = getLogPValueTscoreSpikeTimeAve(spikeDataSet, params, noise_factor_list(nFF));
% %             cellTypeMat(:, nFold, nFF, nTau) = unitGroup;
% %         end
% %     end
% % end
% % cellTypeMatGP43 = cellTypeMat;
% % 
% % % short Ca slow virus
% % nData  = 4;
% % params = DataSetList(nData).params;
% % std_r                  = 0.0246;
% % median_r               = 0.0505;
% % std_d                  = 0.4588;
% % median_d               = 1.7064;
% % load([TempDatDir DataSetList(nData).name '_withOLRemoval.mat'])
% % per_list               = 0.05:0.05:0.95;
% % tau_d_list             = icdf('Normal', per_list, 0, 1) * std_d + median_d;
% % noise_factor_list      = sqrt(icdf('Exponential', per_list, 0.3527));
% % 
% % numNeuron              = length(nDataSet);
% % 
% % cellTypeMat            = nan(numNeuron, numFold, length(noise_factor_list), length(tau_d_list));
% % 
% % for nTau              = 1:length(tau_d_list)
% %     spikeDataSet      = getSyntheticSpikeDeconvDataSimpleVersion(nDataSet, median_r, tau_d_list(nTau), params);  
% %     for nFF           = 1:length(noise_factor_list) 
% %         for nFold     = 1:numFold
% %             unitGroup = getLogPValueTscoreSpikeTimeAve(spikeDataSet, params, noise_factor_list(nFF));
% %             cellTypeMat(:, nFold, nFF, nTau) = unitGroup;
% %         end
% %     end
% % end
% % cellTypeMat6sAAV = cellTypeMat;

load('cellTypeMat.mat')
cellTypeMat = cellTypeMatGP43;
% result for single neuron
figure
nNeuron = 100;
imagesc(squeeze(mean(cellTypeMat(nNeuron, :, :, :), 2)))
axis xy
xlabel('Decay time (sec)')
ylabel('Noise level')

% result for fraction neuron
figure
subplot(1, 2, 1)
imagesc(squeeze(mean(mean(cellTypeMat==1, 1), 2)), [0 1])
colormap('jet')
axis xy
xlabel('Decay time (sec)')
ylabel('Noise level')
title('Frac. Mono. Cell')

subplot(1, 2, 2)
imagesc(squeeze(mean(mean(cellTypeMat==2, 1), 2)), [0 1])
colormap('jet')
axis xy
xlabel('Decay time (sec)')
ylabel('Noise level')
title('Frac. Multi. Cell')

% mean_value
mean(cellTypeMat(:)==1)
mean(cellTypeMat(:)==2)



% result for single neuron
figure
nNeuron = 100;
imagesc(tau_d_list, noise_factor_list, squeeze(mean(cellTypeMat(nNeuron, :, :, :), 2)))
axis xy
xlabel('Decay time (sec)')
ylabel('Noise level')

% result for fraction neuron
figure
subplot(1, 2, 1)
imagesc(tau_d_list, noise_factor_list, squeeze(mean(mean(cellTypeMat==1, 1), 2)), [0 1])
colormap('jet')
axis xy
xlabel('Decay time (sec)')
ylabel('Noise level')
title('Frac. Mono. Cell')

subplot(1, 2, 2)
imagesc(tau_d_list, noise_factor_list, squeeze(mean(mean(cellTypeMat==2, 1), 2)), [0 1])
colormap('jet')
axis xy
xlabel('Decay time (sec)')
ylabel('Noise level')
title('Frac. Multi. Cell')

% mean_value
mean(cellTypeMat(:)==1)

save('cellTypeMat', 'cellTypeMat6sAAV', 'cellTypeMatGP43')
mean(cellTypeMat(:)==2)