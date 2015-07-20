%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 0.  Generate raw activity of all neurons sorted in different ways
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison_V1_0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Common parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('../Func');
setDir;

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
params.Fm              =  20.0000;
params.Kd              =  11.5885;
params.n               =  2.2152;
params.tau_decay       =  1.5471;
params.tau_rise        =  0.0670;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spike
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nDataSet               = getSpikeData(SpikingDataDir, SpikeFileList, params.minNumTrialToAnalysis, params.timeSeries, params.binsize); %#ok<NASGU>
nDataSet               = getDFFSpike(nDataSet, params);
DataSetList(1).name    = 'Shuffle_Spikes';
DataSetList(1).params  = params; 
DataSetList(1).ActiveNeuronIndex = ~findNonActiveNeurons(nDataSet, params);
save([TempDatDir DataSetList(1).name '.mat'], 'nDataSet');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% short Ca fast
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% params.frameRate       =  29.68/2;
% params.binsize         =  1/params.frameRate;
% params.polein          =  -2.6;
% params.poleout         =  -1.4;
% minTimeToAnalysis      =  round(-3.1 * params.frameRate);
% maxTimeToAnalysis      =  round(2.0 * params.frameRate);
% params.timeWindowIndexRange  = minTimeToAnalysis : maxTimeToAnalysis;
% params.timeSeries      = params.timeWindowIndexRange * params.binsize;
% params.minNumTrialToAnalysis =  minNumTrialToAnalysis;
% params.expression      = 'Transgentic';
% nDataSet               = getCaImagingData(CaImagingShortDelayFastDir, CaImagingShortDelayFastFileList, params.minNumTrialToAnalysis, params); %#ok<NASGU>
% nonActiveNeuronIndex   = findNonActiveNeurons(nDataSet, params);
% DataSetList(2).name    = 'Shuffle_Ca_Fast_Short_Delay';
% DataSetList(2).params  = params; 
% DataSetList(2).ActiveNeuronIndex = ~nonActiveNeuronIndex;
% save([TempDatDir DataSetList(2).name '.mat'], 'nDataSet');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% short Ca slow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params.frameRate       =  29.68/2;
params.binsize         =  1/params.frameRate;
params.polein          =  -2.6;
params.poleout         =  -1.4;
minTimeToAnalysis      =  round(-3.1 * params.frameRate);
maxTimeToAnalysis      =  round(2.0 * params.frameRate);
params.timeWindowIndexRange  = minTimeToAnalysis : maxTimeToAnalysis;
params.timeSeries      = params.timeWindowIndexRange * params.binsize;
params.minNumTrialToAnalysis =  minNumTrialToAnalysis;
params.expression      = 'Transgentic';
nDataSet               = getCaImagingData(CaImagingShortDelaySlowDir, CaImagingShortDelaySlowFileList, params.minNumTrialToAnalysis, params); %#ok<NASGU>
nonActiveNeuronIndex   = findNonActiveNeurons(nDataSet, params);
DataSetList(3).name    = 'Shuffle_Ca_Slow_Short_Delay';
DataSetList(3).params  = params; 
DataSetList(3).ActiveNeuronIndex = ~nonActiveNeuronIndex;
save([TempDatDir DataSetList(3).name '.mat'], 'nDataSet');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% long Ca fast
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% params.frameRate       =  29.68/2;
% params.binsize         =  1/params.frameRate;
% params.polein          =  -4.2;
% params.poleout         =  -3.0;
% minTimeToAnalysis      =  round(-4.7 * params.frameRate);
% maxTimeToAnalysis      =  round(2.0 * params.frameRate);
% params.timeWindowIndexRange  = minTimeToAnalysis : maxTimeToAnalysis;
% params.timeSeries      = params.timeWindowIndexRange * params.binsize;
% params.minNumTrialToAnalysis =  minNumTrialToAnalysis;
% params.expression      = 'Transgentic';
% nDataSet               = getCaImagingData(CaImagingLongDelayFastDir, CaImagingLongDelayFastFileList, params.minNumTrialToAnalysis, params); %#ok<NASGU>
% nonActiveNeuronIndex   = findNonActiveNeurons(nDataSet, params);
% DataSetList(4).name    = 'Shuffle_Ca_Fast_Long_Delay';
% DataSetList(4).params  = params; 
% DataSetList(4).ActiveNeuronIndex = ~nonActiveNeuronIndex;
% save([TempDatDir DataSetList(4).name '.mat'], 'nDataSet');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% long Ca slow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% params.frameRate       =  29.68/2;
% params.binsize         =  1/params.frameRate;
% params.polein          =  -4.2;
% params.poleout         =  -3.0;
% minTimeToAnalysis      =  round(-4.7 * params.frameRate);
% maxTimeToAnalysis      =  round(2.0 * params.frameRate);
% params.timeWindowIndexRange  = minTimeToAnalysis : maxTimeToAnalysis;
% params.timeSeries      = params.timeWindowIndexRange * params.binsize;
% params.minNumTrialToAnalysis =  minNumTrialToAnalysis;
% params.expression      = 'Transgentic';
% nDataSet               = getCaImagingData(CaImagingLongDelaySlowDir, CaImagingLongDelaySlowFileList, params.minNumTrialToAnalysis, params); %#ok<NASGU>
% nonActiveNeuronIndex   = findNonActiveNeurons(nDataSet, params);
% DataSetList(5).name    = 'Shuffle_Ca_Slow_Long_Delay';
% DataSetList(5).params  = params; 
% DataSetList(5).ActiveNeuronIndex = ~nonActiveNeuronIndex;
% save([TempDatDir DataSetList(5).name '.mat'], 'nDataSet');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% short Ca slow virus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params.frameRate       =  29.68/2;
params.binsize         =  1/params.frameRate;
params.polein          =  -2.6;
params.poleout         =  -1.4;
minTimeToAnalysis      =  round(-3.1 * params.frameRate);
maxTimeToAnalysis      =  round(2.0 * params.frameRate);
params.timeWindowIndexRange  = minTimeToAnalysis : maxTimeToAnalysis;
params.timeSeries      = params.timeWindowIndexRange * params.binsize;
params.minNumTrialToAnalysis =  minNumTrialToAnalysis;
params.expression      = 'Virus';
nDataSet               = getCaImagingData(CaImagingShortDelaySlowVirusDir, CaImagingShortDelaySlowVirusFileList, params.minNumTrialToAnalysis, params);
nonActiveNeuronIndex   = findNonActiveNeurons(nDataSet, params);
DataSetList(6).name    = 'Shuffle_Ca_Slow_Short_Delay_Virus';
DataSetList(6).params  = params; 
DataSetList(6).ActiveNeuronIndex = ~nonActiveNeuronIndex;
save([TempDatDir DataSetList(6).name '.mat'], 'nDataSet');

save([TempDatDir 'DataListShuffle.mat'], 'DataSetList');

Comparison_V1_0_1;
