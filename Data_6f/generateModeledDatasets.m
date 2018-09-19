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
load([TempDatDir 'Shuffle_Spikes_Nuo_Short_Delay.mat'], 'nDataSet');
params      = DataSetList(1).params;
depth       = [nDataSet.depth_in_um];
spikeDataSet= nDataSet(depth<471);
ActiveNeuronIndex = DataSetList(1).ActiveNeuronIndex;
clear DataSetList

truncatedNormal        = truncate(makedist('Normal'), -0.9, 1.5);
params.Fm              = random(truncatedNormal, length(nDataSet), 1) * 10.6673 + 21.3240;
randPairs              = randi([1 length(nlParaMat)], length(nDataSet), 1);
params.Ca0             = nlParaMat(randPairs, 1);
params.n               = nlParaMat(randPairs, 2);

% Modeled GP5.17
std_r                  = 0.0222;
median_r               = 0.0213;
std_d                  = 0.5390;
median_d               = 0.5898;
params.tau_r           = random(truncatedNormal, length(nDataSet), 1) *  std_r + median_r;
params.tau_d           = random(truncatedNormal, length(nDataSet), 1) *  std_d + median_d;
params.intNoise        = 1.5;
params.extNoise        = 4.5;
nDataSet               = getFakeCaImagingDataSigmoid(spikeDataSet, params);
nData                      = 1;
DataSetList(nData).name    = 'Modeled_GP517';
DataSetList(nData).params  = params; 
DataSetList(nData).ActiveNeuronIndex = ActiveNeuronIndex;
save([TempDatDir DataSetList(nData).name '.mat'], 'nDataSet'); 

save([TempDatDir 'DataListS2C6fModel.mat'], 'DataSetList');