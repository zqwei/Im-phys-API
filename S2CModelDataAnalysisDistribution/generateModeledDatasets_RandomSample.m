%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generating new data using fine tuned params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../Func');
setDir;

% measure of nonlinearity
load([TempDatDir 'FineTunedNLParams.mat'], 'nlParams');

int_noise = [2.0, 2.0];
ext_noise = [0.15, 0.45];

load([TempDatDir 'DataListShuffle.mat'], 'DataSetList');
load([TempDatDir 'Shuffle_Spikes.mat'], 'nDataSet');
params            = DataSetList(1).params;
spikeDataSet      = nDataSet;
ActiveNeuronIndex = DataSetList(1).ActiveNeuronIndex;
clear DataSetList

% Modeled 6s-AAV
nData                  = 1;
nlParamsnData          = squeeze(nlParams(nData, :, :));
truncatedNormal        = truncate(makedist('Normal'), -1.5, 1.5);
% external noise is add after Fm, it does no matter how Fm is given
params.Fm              = nlParamsnData(:, 2);
params.F0              = nlParamsnData(:, 1);
std_r                  = 0.0246;
median_r               = 0.0505;
std_d                  = 0.4588;
median_d               = 1.7064;
params.tau_r           = random(truncatedNormal, length(nDataSet), 1) *  std_r + median_r;
params.intNoise        = int_noise(nData);
params.extNoise        = ext_noise(nData);

mData                  = 0;
Ca0List                = nlParamsnData(:, 3);
taudList               = random(truncatedNormal, length(nDataSet), 1) *  std_d + median_d;

for nParaSet           = 1:1000
    params.Ca0         = prctile(Ca0List, randi(100, length(nDataSet), 1));
    params.n           = 10.^(0.0956-log10(params.Ca0)*0.6010);
    params.tau_d       = prctile(taudList, randi(100, length(nDataSet), 1));

    nDataSet           = getFakeCaImagingDataSigmoidNL(spikeDataSet, params);
    mData              = mData + 1;
    DataSetList(mData).name    = ['Modeled_6s_AAV_nRSParaSet_' num2str(nParaSet, '%04d')];
    DataSetList(mData).params  = params; 
    DataSetList(mData).ActiveNeuronIndex = ActiveNeuronIndex;
    save([TempDatDir DataSetList(mData).name '.mat'], 'nDataSet'); 
end


% Modeled GP43
nData                  = 2;
nlParamsnData          = squeeze(nlParams(nData, :, :));
truncatedNormal        = truncate(makedist('Normal'), -1.5, 1.5);
% external noise is add after Fm, it does no matter how Fm is given
params.Fm              = nlParamsnData(:, 2);
params.F0              = nlParamsnData(:, 1);
std_r                  = 0.0375;
median_r               = 0.0927;
std_d                  = 0.5374;
median_d               = 1.2294;
params.tau_r           = random(truncatedNormal, length(nDataSet), 1) *  std_r + median_r;
params.intNoise        = int_noise(nData);
params.extNoise        = ext_noise(nData);

Ca0List                = nlParamsnData(:, 3);
taudList               = random(truncatedNormal, length(nDataSet), 1) *  std_d + median_d;

for nParaSet           = 1:1000
    params.Ca0         = prctile(Ca0List, randi(100, length(nDataSet), 1));
    params.n           = 10.^(0.0239-log10(params.Ca0)*0.5961);
    params.tau_d       = prctile(taudList, randi(100, length(nDataSet), 1));

    nDataSet           = getFakeCaImagingDataSigmoidNL(spikeDataSet, params);
    mData              = mData + 1;
    DataSetList(mData).name    = ['Modeled_GP43_nRSParaSet_' num2str(nParaSet, '%04d')];
    DataSetList(mData).params  = params; 
    DataSetList(mData).ActiveNeuronIndex = ActiveNeuronIndex;
    save([TempDatDir DataSetList(mData).name '.mat'], 'nDataSet'); 
end

save([TempDatDir 'DataListS2CModelRSDistribution.mat'], 'DataSetList');