%
% TODO
% 06/04/17
% Add code to compute distribution of selectivity for each single neuron
% 
% Here we drop the test of parameter of tau_r
% 
% Comparison based on single unit acitivity
% Generating Ca++ imaging data from ephys data using Tsai-Wen's model
% 
% -------------------------------------------------------------------------
% version 1.0
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Common parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('../Func');
setDir;
load ([TempDatDir 'DataListShuffle.mat']);


% short Ca slow
nData  = 3;
params = DataSetList(nData).params;
std_r                  = 0.0375; %0.0375;
median_r               = 0.0927;
std_d                  = 0.5374; %0.5374;
median_d               = 1.2294;
load([TempDatDir DataSetList(nData).name '_withOLRemoval.mat'])
per_list               = 0.05:0.05:0.95;
tau_d_list             = icdf('Normal', per_list, 0, 1) * std_d + median_d;
noise_factor_list      = sqrt(icdf('Exponential', per_list, 0.3527));

numNeuron              = length(nDataSet);

cellTypeMat            = nan(numNeuron, length(noise_factor_list), length(tau_d_list));

for nTau     = 1:length(tau_d_list)
    
    spikeDataSet = getFakeSpikeDeconvDataSimpleVersion(nDataSet, median_r, tau_d_list(nTau), params);  
    
    for nFF  = 1:length(noise_factor_list) 
        unitGroup = getLogPValueTscoreSpikeTimeAve(spikeDataSet, params, noise_factor_list(nFF));
        cellTypeMat(:, nFF, nTau) = unitGroup;
    end
end

% result for single neuron
figure
nNeuron = 100;
imagesc(tau_d_list, noise_factor_list, squeeze(cellTypeMat(nNeuron,:,:)))
axis xy
xlabel('Decay time (sec)')
ylabel('Noise level')

% result for fraction neuron
figure
subplot(1, 2, 1)
imagesc(tau_d_list, noise_factor_list, squeeze(mean(cellTypeMat==1)), [0 1])
colormap('jet')
axis xy
xlabel('Decay time (sec)')
ylabel('Noise level')
title('Frac. Mono. Cell')

subplot(1, 2, 2)
imagesc(tau_d_list, noise_factor_list, squeeze(mean(cellTypeMat==2)), [0 1])
colormap('jet')
axis xy
xlabel('Decay time (sec)')
ylabel('Noise level')
title('Frac. Multi. Cell')

% mean_value
mean(cellTypeMat(:)==1)
mean(cellTypeMat(:)==2)


% short Ca slow virus
nData  = 4;
params = DataSetList(nData).params;
std_r                  = 0.0246;
median_r               = 0.0505;
std_d                  = 0.4588;
median_d               = 1.7064;
load([TempDatDir DataSetList(nData).name '_withOLRemoval.mat'])
per_list               = 0.05:0.05:0.95;
tau_d_list             = icdf('Normal', per_list, 0, 1) * std_d + median_d;
noise_factor_list      = sqrt(icdf('Exponential', per_list, 0.3527));

numNeuron              = length(nDataSet);

cellTypeMat            = nan(numNeuron, length(noise_factor_list), length(tau_d_list));

for nTau     = 1:length(tau_d_list)
    
    spikeDataSet = getFakeSpikeDeconvDataSimpleVersion(nDataSet, median_r, tau_d_list(nTau), params);  
    
    for nFF  = 1:length(noise_factor_list) 
        unitGroup = getLogPValueTscoreSpikeTimeAve(spikeDataSet, params, noise_factor_list(nFF));
        cellTypeMat(:, nFF, nTau) = unitGroup;
    end
end

% result for single neuron
figure
nNeuron = 100;
imagesc(tau_d_list, noise_factor_list, squeeze(cellTypeMat(nNeuron,:,:)))
axis xy
xlabel('Decay time (sec)')
ylabel('Noise level')

% result for fraction neuron
figure
subplot(1, 2, 1)
imagesc(tau_d_list, noise_factor_list, squeeze(mean(cellTypeMat==1)), [0 1])
colormap('jet')
axis xy
xlabel('Decay time (sec)')
ylabel('Noise level')
title('Frac. Mono. Cell')

subplot(1, 2, 2)
imagesc(tau_d_list, noise_factor_list, squeeze(mean(cellTypeMat==2)), [0 1])
colormap('jet')
axis xy
xlabel('Decay time (sec)')
ylabel('Noise level')
title('Frac. Multi. Cell')

% mean_value
mean(cellTypeMat(:)==1)
mean(cellTypeMat(:)==2)
