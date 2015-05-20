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

clear all %#ok<CLSCR>

addpath('../Func');
setDir;
load([TempDatDir 'DataListModeled.mat'], 'DataSetList');
oldDataSetList         = DataSetList;
clear DataSetList;

minNumTrialToAnalysis  = 20;
nData                  = 0;
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

%
% Raw activity
%
% Spike
% nDataSet               = getSpikeDataWithEphysTime(SpikingDataDir, SpikeFileList, params.minNumTrialToAnalysis, params.timeSeries, params.binsize);
% DataSetList(1).name    = 'Modeled_Spikes';
% DataSetList(1).params  = params; 
% save([TempDatDir DataSetList(1).name '.mat'], 'nDataSet');
nData                      = nData + 1;
load([TempDatDir oldDataSetList(nData).name '.mat'], 'nDataSet');
spikeDataSet           = nDataSet;
% nDataSet               = getDFFSpike(spikeDataSet, params);
% nData                      = nData + 1;
DataSetList(nData).name    = 'Modeled_Spikes';
DataSetList(nData).params  = params; 
% save([TempDatDir DataSetList(nData).name '.mat'], 'nDataSet');


% Modeled GCaMP6s Virus
truncatedNormal        = truncate(makedist('Normal'), -1.5, 1.5);
params.Fm              = random(truncatedNormal, length(nDataSet), 1) * 10.6673 + 21.3240;
params.K               = random(truncatedNormal, length(nDataSet), 1) *  5.8337 + 13.9248;
params.n               = random(truncatedNormal, length(nDataSet), 1) *  0.4106 +  1.7531;
params.tau_r           = random(truncatedNormal, length(nDataSet), 1) *  0.0344 +  0.0728;
params.tau_d           = random(truncatedNormal, length(nDataSet), 1) *  0.5055 +  1.4551;
params.tau_r           = params.tau_r;
params.tau_d           = params.tau_d;
params.intNoise        = 0.0;
params.extNoise        = 0.0;
nDataSet               = getFakeCaImagingData(spikeDataSet, params);
nData                      = nData + 1;
DataSetList(nData).name    = 'Modeled_Ca_Long_Decay_No_Noise';
DataSetList(nData).params  = params; 
save([TempDatDir DataSetList(nData).name '.mat'], 'nDataSet'); 



% Modeled GCaMP6s Virus
% truncatedNormal        = truncate(makedist('Normal'), -1.5, 1.5);
% params.Fm              = random(truncatedNormal, length(nDataSet), 1) * 10.6673 + 21.3240;
% params.K               = random(truncatedNormal, length(nDataSet), 1) *  5.8337 + 13.9248;
% params.n               = random(truncatedNormal, length(nDataSet), 1) *  0.4106 +  1.7531;
% params.tau_r           = random(truncatedNormal, length(nDataSet), 1) *  0.0344 +  0.0728;
% params.tau_d           = random(truncatedNormal, length(nDataSet), 1) *  0.5055 +  1.4551;
% params.tau_r           = params.tau_r;
% params.tau_d           = params.tau_d/5;
% params.intNoise        = 0.0;
% params.extNoise        = 0.0;
% nDataSet               = getFakeCaImagingData(spikeDataSet, params);
% nData                      = nData + 1;
% DataSetList(nData).name    = 'Modeled_Ca_Short_Decay_No_Noise';
% DataSetList(nData).params  = params; 
% save([TempDatDir DataSetList(nData).name '.mat'], 'nDataSet'); 


% % Modeled GCaMP6s Virus
% truncatedNormal        = truncate(makedist('Normal'), -1.5, 1.5);
% params.Fm              = random(truncatedNormal, length(nDataSet), 1) * 10.6673 + 21.3240;
% params.K               = random(truncatedNormal, length(nDataSet), 1) *  5.8337 + 13.9248;
% params.n               = random(truncatedNormal, length(nDataSet), 1) *  0.4106 +  1.7531;
% params.tau_r           = random(truncatedNormal, length(nDataSet), 1) *  0.0344 +  0.0728;
% params.tau_d           = random(truncatedNormal, length(nDataSet), 1) *  0.5055 +  1.4551;
% params.tau_r           = params.tau_r;
% params.tau_d           = params.tau_r;
% params.intNoise        = 0.0;
% params.extNoise        = 0.0;
% nDataSet               = getFakeCaImagingData(spikeDataSet, params);
% DataSetList(4).name    = 'Modeled_Ca_VeryShort_Decay_No_Noise';
% DataSetList(4).params  = params; 
% save([TempDatDir DataSetList(4).name '.mat'], 'nDataSet'); 


% Modeled GCaMP6s Virus
truncatedNormal        = truncate(makedist('Normal'), -1.5, 1.5);
params.Fm              = random(truncatedNormal, length(nDataSet), 1) * 10.6673 + 21.3240;
params.K               = random(truncatedNormal, length(nDataSet), 1) *  5.8337 + 13.9248;
params.n               = random(truncatedNormal, length(nDataSet), 1) *  0.4106 +  1.7531;
params.tau_r           = random(truncatedNormal, length(nDataSet), 1) *  0.0344 +  0.0728;
params.tau_d           = random(truncatedNormal, length(nDataSet), 1) *  0.5055 +  1.4551;
params.tau_r           = params.tau_r;
params.tau_d           = params.tau_d;
params.intNoise        = 0.5;
params.extNoise        = 0.0;
nDataSet               = getFakeCaImagingData(spikeDataSet, params);
nData                      = nData + 1;
DataSetList(nData).name    = 'Modeled_Ca_Long_Decay_Const_Noise';
DataSetList(nData).params  = params; 
save([TempDatDir DataSetList(nData).name '.mat'], 'nDataSet'); 
% 
% 
% % Modeled GCaMP6s Virus
% truncatedNormal        = truncate(makedist('Normal'), -1.5, 1.5);
% params.Fm              = random(truncatedNormal, length(nDataSet), 1) * 10.6673 + 21.3240;
% params.K               = random(truncatedNormal, length(nDataSet), 1) *  5.8337 + 13.9248;
% params.n               = random(truncatedNormal, length(nDataSet), 1) *  0.4106 +  1.7531;
% params.tau_r           = random(truncatedNormal, length(nDataSet), 1) *  0.0344 +  0.0728;
% params.tau_d           = random(truncatedNormal, length(nDataSet), 1) *  0.5055 +  1.4551;
% params.tau_r           = params.tau_r;
% params.tau_d           = params.tau_d/5;
% params.linearNoise     = 0.4;
% params.constNoise      = 0.0;
% nDataSet               = getFakeCaImagingData(spikeDataSet, params);
% DataSetList(5).name    = 'Modeled_Ca_Short_Decay_Linear_Noise';
% DataSetList(5).params  = params; 
% save([TempDatDir DataSetList(5).name '.mat'], 'nDataSet'); 
% 
% 
% % Modeled GCaMP6s Virus
% truncatedNormal        = truncate(makedist('Normal'), -1.5, 1.5);
% params.Fm              = random(truncatedNormal, length(nDataSet), 1) * 10.6673 + 21.3240;
% params.K               = random(truncatedNormal, length(nDataSet), 1) *  5.8337 + 13.9248;
% params.n               = random(truncatedNormal, length(nDataSet), 1) *  0.4106 +  1.7531;
% params.tau_r           = random(truncatedNormal, length(nDataSet), 1) *  0.0344 +  0.0728;
% params.tau_d           = random(truncatedNormal, length(nDataSet), 1) *  0.5055 +  1.4551;
% params.tau_r           = params.tau_r;
% params.tau_d           = params.tau_d/10;
% params.linearNoise     = 0.0;
% params.constNoise      = 0.4;
% nDataSet               = getFakeCaImagingData(spikeDataSet, params);
% DataSetList(6).name    = 'Modeled_Ca_VeryShort_Decay_Const_Noise';
% DataSetList(6).params  = params; 
% save([TempDatDir DataSetList(6).name '.mat'], 'nDataSet'); 

save([TempDatDir 'DataListModeled.mat'], 'DataSetList');

Comparison_V1_0_1;
Comparison_V1_0_2;