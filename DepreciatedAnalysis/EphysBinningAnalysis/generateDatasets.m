%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 0.  Generate raw activity of all neurons sorted in different ways
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison_V2_0

addpath('../Func');
setDir;

minNumTrialToAnalysis  = 20;

% Shuffle_Spikes_non_overlap
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
minFiringRate          = 5; % Hz per epoch
nDataSet               = getSpikeDataWithEphysTime(SpikingDataDir, SpikeFileList, params.minNumTrialToAnalysis, params.timeSeries, params.binsize);                                  
DataSetList(1).name    = 'Shuffle_Spikes_non_overlap';
DataSetList(1).params  = params; 
DataSetList(1).ActiveNeuronIndex = findHighFiringUnits(nDataSet, params, minFiringRate);
save([TempDatDir DataSetList(1).name '.mat'], 'nDataSet');


params.binsize         =  0.070;
params.frameRate       =  200; %Hz -- smallest bin is 0.001
minTimeToAnalysis      =  round(-3.1 * params.frameRate);
maxTimeToAnalysis      =  round(2.0 * params.frameRate);
params.timeWindowIndexRange  = minTimeToAnalysis : maxTimeToAnalysis;
params.timeSeries      = params.timeWindowIndexRange / params.frameRate;
nDataSet               =  getSpikeDataWithEphysTimeAndBoxCar(SpikingDataDir, SpikeFileList, params.minNumTrialToAnalysis, params.timeSeries, params.binsize);                                  
DataSetList(2).name    = 'Shuffle_Spikes_boxcar_070ms';
DataSetList(2).params  = params; 
DataSetList(2).ActiveNeuronIndex = findHighFiringUnits(nDataSet, params, minFiringRate);
save([TempDatDir DataSetList(2).name '.mat'], 'nDataSet', '-v7.3');


params.binsize         =  0.100;
nDataSet               =  getSpikeDataWithEphysTimeAndBoxCar(SpikingDataDir, SpikeFileList, params.minNumTrialToAnalysis, params.timeSeries, params.binsize);                                  
DataSetList(3).name    = 'Shuffle_Spikes_boxcar_100ms';
DataSetList(3).params  = params; 
DataSetList(3).ActiveNeuronIndex = findHighFiringUnits(nDataSet, params, minFiringRate);
save([TempDatDir DataSetList(3).name '.mat'], 'nDataSet', '-v7.3');

params.binsize         =  0.150;
nDataSet               =  getSpikeDataWithEphysTimeAndBoxCar(SpikingDataDir, SpikeFileList, params.minNumTrialToAnalysis, params.timeSeries, params.binsize);                                  
DataSetList(4).name    = 'Shuffle_Spikes_boxcar_150ms';
DataSetList(4).params  = params; 
DataSetList(4).ActiveNeuronIndex = findHighFiringUnits(nDataSet, params, minFiringRate);
save([TempDatDir DataSetList(4).name '.mat'], 'nDataSet', '-v7.3');

params.binsize         =  0.200;
nDataSet               =  getSpikeDataWithEphysTimeAndBoxCar(SpikingDataDir, SpikeFileList, params.minNumTrialToAnalysis, params.timeSeries, params.binsize);                                  
DataSetList(5).name    = 'Shuffle_Spikes_boxcar_200ms';
DataSetList(5).params  = params; 
DataSetList(5).ActiveNeuronIndex = findHighFiringUnits(nDataSet, params, minFiringRate);
save([TempDatDir DataSetList(5).name '.mat'], 'nDataSet', '-v7.3');

params.binsize         =  0.250;
nDataSet               =  getSpikeDataWithEphysTimeAndBoxCar(SpikingDataDir, SpikeFileList, params.minNumTrialToAnalysis, params.timeSeries, params.binsize);                                  
DataSetList(6).name    = 'Shuffle_Spikes_boxcar_250ms';
DataSetList(6).params  = params; 
DataSetList(6).ActiveNeuronIndex = findHighFiringUnits(nDataSet, params, minFiringRate);
save([TempDatDir DataSetList(6).name '.mat'], 'nDataSet', '-v7.3');



save([TempDatDir 'DataListEphys.mat'], 'DataSetList');
