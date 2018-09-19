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
load([TempDatDir 'nlParaMat.mat'],'nlParaMat');
clear params;
load([TempDatDir 'DataListShuffle.mat'], 'DataSetList');
load([TempDatDir 'Shuffle_Spikes.mat'], 'nDataSet');
params      = DataSetList(1).params;
spikeDataSet           = nDataSet;
ActiveNeuronIndex = DataSetList(1).ActiveNeuronIndex;
clear DataSetList

% Modeled 6s-AAV
truncatedNormal        = truncate(makedist('Normal'), -1.5, 1.5);
% external noise is add after Fm, it does no matter how Fm is given
params.Fm              = random(truncatedNormal, length(nDataSet), 1) * 10.6673 + 21.3240;
% std_K                  = min(S2Cparams(idK).K/3, S2Cparams(12).K);
% Ca0 and n are given as pairs (randomly sampled)
randPairs              = randi([1 length(nlParaMat)], length(nDataSet), 1);
params.Ca0             = nlParaMat(randPairs, 1);
params.n               = nlParaMat(randPairs, 2);
% tau_r and tau_d are given randomly
std_r                  = 0.0246;
median_r               = 0.0505;
std_d                  = 0.4588;
median_d               = 1.7064;
params.tau_r           = random(truncatedNormal, length(nDataSet), 1) *  std_r + median_r;
params.tau_d           = random(truncatedNormal, length(nDataSet), 1) *  std_d + median_d;
params.intNoise        = 1.5;
params.extNoise        = 1.5;
nDataSet               = getFakeCaImagingDataSigmoid(spikeDataSet, params);
nData                      = 1;
DataSetList(nData).name    = 'Modeled_6s_AAV';
DataSetList(nData).params  = params; 
DataSetList(nData).ActiveNeuronIndex = ActiveNeuronIndex;
save([TempDatDir DataSetList(nData).name '.mat'], 'nDataSet'); 

% Modeled GP4.3
% tau_r and tau_d are given randomly
std_r                  = 0.0375;
median_r               = 0.0927;
std_d                  = 0.5374;
median_d               = 1.2294;
params.tau_r           = random(truncatedNormal, length(nDataSet), 1) *  std_r + median_r;
params.tau_d           = random(truncatedNormal, length(nDataSet), 1) *  std_d + median_d;
params.intNoise        = 1.5;
params.extNoise        = 4.5;
nDataSet               = getFakeCaImagingDataSigmoid(spikeDataSet, params);
nData                      = 2;
DataSetList(nData).name    = 'Modeled_GP43';
DataSetList(nData).params  = params; 
DataSetList(nData).ActiveNeuronIndex = ActiveNeuronIndex;
save([TempDatDir DataSetList(nData).name '.mat'], 'nDataSet'); 

save([TempDatDir 'DataListS2CModel.mat'], 'DataSetList');