%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Common parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The data is F not dF/F 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('../Func');
setDir;

minNumTrialToAnalysis  = 20;

% Set data directory
minRate                = 5;
perMinRate             = 0.4;
ROCThres               = 0.7;
minUnitsSession        = 3;

% short Ca fast
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
nDataSet               = getCaImagingDataRaw(CaImagingShortDelayFastDir, CaImagingShortDelayFastFileList, params.minNumTrialToAnalysis, params); %#ok<NASGU>
DataSetList(2).name    = 'Raw_Ca_Fast_Short_Delay';
DataSetList(2).params  = params; 
save([TempDatDir DataSetList(2).name '.mat'], 'nDataSet');

% short Ca slow
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
nDataSet               = getCaImagingDataRaw(CaImagingShortDelaySlowDir, CaImagingShortDelaySlowFileList, params.minNumTrialToAnalysis, params); %#ok<NASGU>
DataSetList(3).name    = 'Raw_Ca_Slow_Short_Delay';
DataSetList(3).params  = params; 
save([TempDatDir DataSetList(3).name '.mat'], 'nDataSet');

% short Ca slow
params.frameRate       =  29.68/2;
params.binsize         =  1/params.frameRate;
params.polein          =  -4.2;
params.poleout         =  -3.0;
minTimeToAnalysis      =  round(-4.7 * params.frameRate);
maxTimeToAnalysis      =  round(2.0 * params.frameRate);
params.timeWindowIndexRange  = minTimeToAnalysis : maxTimeToAnalysis;
params.timeSeries      = params.timeWindowIndexRange * params.binsize;
params.minNumTrialToAnalysis =  minNumTrialToAnalysis;
params.expression      = 'Transgentic';
nDataSet               = getCaImagingDataRaw(CaImagingLongDelayFastDir, CaImagingLongDelayFastFileList, params.minNumTrialToAnalysis, params); %#ok<NASGU>
DataSetList(4).name    = 'Raw_Ca_Fast_Long_Delay';
DataSetList(4).params  = params; 
save([TempDatDir DataSetList(4).name '.mat'], 'nDataSet');

% long Ca slow
params.frameRate       =  29.68/2;
params.binsize         =  1/params.frameRate;
params.polein          =  -4.2;
params.poleout         =  -3.0;
minTimeToAnalysis      =  round(-4.7 * params.frameRate);
maxTimeToAnalysis      =  round(2.0 * params.frameRate);
params.timeWindowIndexRange  = minTimeToAnalysis : maxTimeToAnalysis;
params.timeSeries      = params.timeWindowIndexRange * params.binsize;
params.minNumTrialToAnalysis =  minNumTrialToAnalysis;
params.expression      = 'Transgentic';
nDataSet               = getCaImagingDataRaw(CaImagingLongDelaySlowDir, CaImagingLongDelaySlowFileList, params.minNumTrialToAnalysis, params); %#ok<NASGU>
DataSetList(5).name    = 'Raw_Ca_Slow_Long_Delay';
DataSetList(5).params  = params; 
save([TempDatDir DataSetList(5).name '.mat'], 'nDataSet');


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
nDataSet               = getCaImagingDataRaw(CaImagingShortDelaySlowVirusDir, CaImagingShortDelaySlowVirusFileList, params.minNumTrialToAnalysis, params);
DataSetList(6).name    = 'Raw_Ca_Slow_Short_Delay_Virus';
DataSetList(6).params  = params; 
save([TempDatDir DataSetList(6).name '.mat'], 'nDataSet');

save([TempDatDir 'DataListShuffleRaw.mat'], 'DataSetList');
