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
load('ParamsFitCells_S2CModel_Sim.mat');
S2Cparams   = params;
clear params;
idTauD      = 1;
idN         = 1;
idK         = 1;
minNumTrialToAnalysis  = 20;
params.frameRate       =  29.68/2;
params.binsize         =  1/params.frameRate;
params.polein          =  -2.6;
params.poleout         =  -1.3;
minTimeToAnalysis      =  round(-3.1 * params.frameRate);
maxTimeToAnalysis      =  round(2.0 * params.frameRate);
params.timeWindowIndexRange  = minTimeToAnalysis : maxTimeToAnalysis;
params.timeSeries      = params.timeWindowIndexRange * params.binsize;
params.minNumTrialToAnalysis =  minNumTrialToAnalysis;
params.expression      = 'None';
minFiringRate          = 5;
% Raw activity
% Spike
nDataSet               = getSpikeDataWithEphysTime(SpikingDataDir, SpikeFileList, params.minNumTrialToAnalysis, params.timeSeries, params.binsize);
spikeDataSet           = nDataSet;
ActiveNeuronIndex = findHighFiringUnits(spikeDataSet, params, minFiringRate);


% Modeled GCaMP6s long decay
truncatedNormal        = truncate(makedist('Normal'), -1.5, 1.5);
params.Fm              = random(truncatedNormal, length(nDataSet), 1) * 10.6673 + 21.3240;
std_K                  = min(S2Cparams(idK).K/3, S2Cparams(12).K);
params.K               = random(truncatedNormal, length(nDataSet), 1) *  std_K               + S2Cparams(idN).K;
params.n               = random(truncatedNormal, length(nDataSet), 1) *  S2Cparams(12).n     + S2Cparams(idN).n;
params.tau_r           = random(truncatedNormal, length(nDataSet), 1) *  S2Cparams(12).tau_r + S2Cparams(1).tau_r;
params.tau_d           = random(truncatedNormal, length(nDataSet), 1) *  S2Cparams(12).tau_d + S2Cparams(idTauD).tau_d;
params.Fm              = params.Fm;
params.tau_r           = params.tau_r;
params.tau_d           = params.tau_d;
params.intNoise        = 1.5;
params.extNoise        = 1.5;
nDataSet               = getFakeCaImagingData(spikeDataSet, params);
nData                      = 1;
DataSetList(nData).name    = 'Modeled_Ca_Long_Decay_No_Noise';
DataSetList(nData).params  = params; 
DataSetList(nData).ActiveNeuronIndex = ActiveNeuronIndex;
save([TempDatDir DataSetList(nData).name '.mat'], 'nDataSet'); 

% Modeled GCaMP6s short decay
params.Fm              = params.Fm;
params.tau_r           = params.tau_r;
params.tau_d           = params.tau_d*0.5;
params.intNoise        = 1.5;
params.extNoise        = 1.5;
nDataSet               = getFakeCaImagingData(spikeDataSet, params);
nData                      = 2;
DataSetList(nData).name    = 'Modeled_Ca_Short_Decay_No_Noise';
DataSetList(nData).params  = params; 
DataSetList(nData).ActiveNeuronIndex = ActiveNeuronIndex;
save([TempDatDir DataSetList(nData).name '.mat'], 'nDataSet'); 

save([TempDatDir 'DataListS2CModel.mat'], 'DataSetList');