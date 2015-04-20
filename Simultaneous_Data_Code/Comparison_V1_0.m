% 
% Comparison based on single unit acitivity
% 
% -------------------------------------------------------------------------
% version 1.0
%
% Comparison list
%
% 1.  Raw activity of all neurons sorted in different ways
% 2.  Shape of activity distribution across trial periods (pre-sample,
%     sample, delay, response)
% 3.  Selectivity over time
% 4.  Single neuron choice probability index distribution (ROC), across
%     trial periods
% 5.  Per session decodability over time
% 6.  Collected population decision decodability over time
% 7.  Saturation of decodeability with number of neurons
% 8.  Decoding which trial period one is in, how fast can you tell the data
%     that the trial period switched
% 9.  Measures of trial-to-trial variability, by trial period
% 10. Fraction of variance captured by PCA, per session, by trial period.
%
% -------------------------------------------------------------------------
% version 1.1
% 
% + Save the dataset while before rerun all the code. 
% Simultaneously recording data
% 
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1.  Raw activity of all neurons sorted in different ways
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison_V1_1
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.  Shape of activity distribution across trial periods (pre-sample,
%     sample, delay, response)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison_V1_2
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3.  Selectivity over time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison_V1_3
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4.  Single neuron choice probability index distribution (ROC), across
%     trial periods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison_V1_4
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5.  Per session decodability over time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison_V1_5
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 6.  Collected population decision decodability over time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison_V1_6
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 7.  Saturation of decodeability with number of neurons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison_V1_7
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 8.  Decoding which trial period one is in, how fast can you tell the data
%     that the trial period switched
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison_V1_8
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 9.  Measures of trial-to-trial variability, by trial period
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison_V1_9
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 10. Fraction of variance captured by PCA, per session, by trial period.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison_V1_10
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Common parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('../Func');
setDir;

params.frameRate       =  29.68/2;
params.binsize         =  1/params.frameRate;
params.polein          =  -2.6;
params.poleout         =  -1.3;
minTimeToAnalysis      =  round(-3.1 * params.frameRate);
maxTimeToAnalysis      =  round(2.0 * params.frameRate);
params.timeWindowIndexRange  = minTimeToAnalysis : maxTimeToAnalysis;
params.timeSeries      = params.timeWindowIndexRange * params.binsize;
params.minNumTrialToAnalysis =  30;
params.Fm              =  20.0000;
params.Kd              =  11.5885;
params.n               =  2.2152;
params.tau_decay       =  1.5471;
params.tau_rise        =  0.0670;


minRate                = 5;
perMinRate             = 0.4;
ROCThres               = 0.70;
minUnitsSession        = 3;
%
% Raw activity

% Spike
nDataSet               = getSpikeData(SpikingDataDir, SpikeFileList, params.minNumTrialToAnalysis, params.timeSeries, params.binsize);
SpikeDataSet           = nDataSet;
[nDataSet3D, nDataSet] = getSimultaneousSpikeData(nDataSet, params, minRate, perMinRate, ROCThres, minUnitsSession); %#ok<ASGLU>
DataSetList(1).name    = 'Spikes';
DataSetList(1).params  = params; 
save([TempDatDir DataSetList(1).name '.mat'], 'nDataSet', 'nDataSet3D');



minRate                = 5;
perMinRate             = 0.4;
ROCThres               = 0.7;
minUnitsSession        = 3;

% short Ca fast
params.poleout         =  -1.4;
nDataSet               = getCaImagingData(CaImagingShortDelayFastDir, CaImagingShortDelayFastFileList, params.minNumTrialToAnalysis, params);
[nDataSet3D, nDataSet] = getSimultaneousCaimagingData(nDataSet, params, ROCThres, minUnitsSession); %#ok<ASGLU>
DataSetList(2).name    = 'Ca_Fast_Short_Delay';
DataSetList(2).params  = params; 
save([TempDatDir DataSetList(2).name '.mat'], 'nDataSet', 'nDataSet3D');

% short Ca slow
params.poleout         =  -1.4;
nDataSet               = getCaImagingData(CaImagingShortDelaySlowDir, CaImagingShortDelaySlowFileList, params.minNumTrialToAnalysis, params);
[nDataSet3D, nDataSet] = getSimultaneousCaimagingData(nDataSet, params, ROCThres, minUnitsSession); %#ok<ASGLU>
DataSetList(3).name    = 'Ca_Slow_Short_Delay';
DataSetList(3).params  = params; 
save([TempDatDir DataSetList(3).name '.mat'], 'nDataSet', 'nDataSet3D');

% long Ca fast
params.polein          =  -4.2;
params.poleout         =  -3.0;
minTimeToAnalysis      =  round(-4.7 * params.frameRate);
maxTimeToAnalysis      =  round(2.0 * params.frameRate);
params.timeWindowIndexRange  = minTimeToAnalysis : maxTimeToAnalysis;
params.timeSeries      = params.timeWindowIndexRange * params.binsize;
nDataSet               = getCaImagingData(CaImagingLongDelayFastDir, CaImagingLongDelayFastFileList, params.minNumTrialToAnalysis, params);
[nDataSet3D, nDataSet] = getSimultaneousCaimagingData(nDataSet, params, ROCThres, minUnitsSession); %#ok<ASGLU>
DataSetList(4).name    = 'Ca_Fast_Long_Delay';
DataSetList(4).params  = params; 
save([TempDatDir DataSetList(4).name '.mat'], 'nDataSet', 'nDataSet3D');

% long Ca slow
nDataSet               = getCaImagingData(CaImagingLongDelaySlowDir, CaImagingLongDelaySlowFileList, params.minNumTrialToAnalysis, params);
[nDataSet3D, nDataSet] = getSimultaneousCaimagingData(nDataSet, params, ROCThres, minUnitsSession);
DataSetList(5).name    = 'Ca_Slow_Long_Delay';
DataSetList(5).params  = params; 
save([TempDatDir DataSetList(5).name '.mat'], 'nDataSet', 'nDataSet3D');

% % Fake Ca
% nDataSet               = getFakeCaImagingData(SpikingDataDir, SpikeFileList, params.minNumTrialToAnalysis, params);
% [nDataSet3D, nDataSet] = getSimultaneousFakeCaData(nDataSet, SpikeDataSet, params, minRate, perMinRate, ROCThres, minUnitsSession); %#ok<ASGLU>
% DataSetList(2).name    = 'Fake_Ca';
% DataSetList(2).params  = params; 
% save([TempDatDir DataSetList(2).name '.mat'], 'nDataSet', 'nDataSet3D');
% 

save([TempDatDir 'DataList.mat'], 'DataSetList');

