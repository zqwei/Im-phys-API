%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Selectiviy index of Nuo's version
% Nomarilized contra - ispi firing rate
% 
% For this analysis, firing rate data for spiking neurons are recomputed
% using Nuo's filter
% 
% Those for Ca++ remain the same as the other analysis
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../Func');
setDir;
load ([TempDatDir 'DataListModeled.mat']);

%%%% For spiking data

if ~exist([PlotDir 'ModeledSingleUnitsNormalizedDiff'],'dir')
    mkdir([PlotDir 'ModeledSingleUnitsNormalizedDiff'])
end

params.frameRate       =  1000;
params.binsize         =  1/params.frameRate; % 1 ms bin
params.polein          =  -2.6;
params.poleout         =  -1.3;
minTimeToAnalysis      =  round(-3.5 * params.frameRate); % these parameters from Nuo
maxTimeToAnalysis      =  round(2.0 * params.frameRate); % these parameters from Nuo
params.timeWindowIndexRange  = minTimeToAnalysis : maxTimeToAnalysis;
params.timeSeries      = params.timeWindowIndexRange * params.binsize;
params.minNumTrialToAnalysis =  20;
params.expression      = 'None';
nDataSet               = getSpikeData(SpikingDataDir, SpikeFileList, params.minNumTrialToAnalysis, params.timeSeries, params.binsize);

cellType               = false(length(nDataSet), 1);

for nUnit  = 1:length(nDataSet)
    if strcmp(nDataSet(nUnit).cell_type, 'putative_interneuron')
        cellType(nUnit)   = false;
    elseif strcmp(nDataSet(nUnit).cell_type, 'putative_pyramidal')
        cellType(nUnit)   = true;
    end   
end

pryCellIndex    = cellType;
nDataSet        = nDataSet(pryCellIndex);
PSTHStartTime   = -3.5;
PSTHEndTime     = 2.2;
wholeTrialStartingTime = -3.1; % pole in = -2.6 % pole out = -1.4
wholeTrialEndTime      = 1.3; % these two values are from Nuo
plotSelectivityIndex(nDataSet, PSTHStartTime, PSTHEndTime, wholeTrialStartingTime, wholeTrialEndTime, params);

setPrint(16, 6*2*2, [PlotDir 'ModeledSingleUnitsNormalizedDiff/SingleUnitsNormalizedDiff_' DataSetList(1).name], 'pdf')

for nData = 2:length(DataSetList)    
    load([TempDatDir DataSetList(nData).name '.mat'])
    plotCaSelectivityIndex(nDataSet(pryCellIndex), wholeTrialStartingTime, wholeTrialEndTime, DataSetList(nData).params);    
    setPrint(16, 6*2*2, [PlotDir 'ModeledSingleUnitsNormalizedDiff/SingleUnitsNormalizedDiff_' DataSetList(nData).name], 'pdf')
end

close all