% 
% Comparison based on single unit acitivity
% 
% -------------------------------------------------------------------------
% version 1.0
%
% Comparison list
%
% 1.  Reproduce Figures in Nuo's paper (Li et al., 2015, Nature)
% 
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Common parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('../Func');
setDir;

params.frameRate       =  1000;
params.binsize         =  1/params.frameRate; % 1 ms bin
params.polein          =  -2.6;
params.poleout         =  -1.3;
minTimeToAnalysis      =  round(-3.5 * params.frameRate); % these parameters from Nuo
maxTimeToAnalysis      =  round(2.0 * params.frameRate); % these parameters from Nuo
params.timeWindowIndexRange  = minTimeToAnalysis : maxTimeToAnalysis;
params.timeSeries      = params.timeWindowIndexRange * params.binsize;
params.minNumTrialToAnalysis =  0;
params.expression      = 'None';


% Spike
nDataSet               = getSpikeData(SpikingDataDir, SpikeFileList, params.minNumTrialToAnalysis, params.timeSeries, params.binsize);
DataSetList(1).name    = 'LiAnalysis_Spikes';
DataSetList(1).params  = params; 
save([TempDatDir DataSetList(1).name '.mat'], 'nDataSet','-v7.3');


save([TempDatDir 'LiAnalysis_DataList.mat'], 'DataSetList');

