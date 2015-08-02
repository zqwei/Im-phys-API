%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 0.  Generate raw activity of all neurons sorted in different ways
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison_V2_0

addpath('../Func');
setDir;

minNumTrialToAnalysis  = 20;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spike
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
nDataSet               = getSpikeData(SpikingDataDir, SpikeFileList, ...
                                      params.minNumTrialToAnalysis, ...
                                      params.timeSeries, params.binsize);
DataSetList(1).name    = 'Shuffle_Spikes';
DataSetList(1).params  = params; 
DataSetList(1).ActiveNeuronIndex = findHighFiringUnits(nDataSet, params, minFiringRate);
save([TempDatDir DataSetList(1).name '.mat'], 'nDataSet');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spike -- DF/F
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nDataSet               = getDFFSpike(nDataSet, params);
DataSetList(7).name    = 'Shuffle_DFFSpikes';
DataSetList(7).params  = params; 
DataSetList(7).ActiveNeuronIndex = ~findNonActiveNeurons(nDataSet, params) & DataSetList(1).ActiveNeuronIndex;
save([TempDatDir DataSetList(7).name '.mat'], 'nDataSet');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% short Ca fast
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
nDataSet               = getCaImagingData(CaImagingShortDelayFastDir, ...
                                          CaImagingShortDelayFastFileList, ...
                                          params.minNumTrialToAnalysis, params); 
nonActiveNeuronIndex   = findNonActiveNeurons(nDataSet, params);
DataSetList(2).name    = 'Shuffle_Ca_Fast_Short_Delay';
DataSetList(2).params  = params; 
DataSetList(2).ActiveNeuronIndex = ~nonActiveNeuronIndex;
save([TempDatDir DataSetList(2).name '.mat'], 'nDataSet');

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
nDataSet               = getCaImagingData(CaImagingShortDelaySlowDir, ...
                                          CaImagingShortDelaySlowFileList, ...
                                          params.minNumTrialToAnalysis, params); 
nonActiveNeuronIndex   = findNonActiveNeurons(nDataSet, params);
DataSetList(3).name    = 'Shuffle_Ca_Slow_Short_Delay';
DataSetList(3).params  = params; 
DataSetList(3).ActiveNeuronIndex = ~nonActiveNeuronIndex;
save([TempDatDir DataSetList(3).name '.mat'], 'nDataSet');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% long Ca fast
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
nDataSet               = getCaImagingData(CaImagingLongDelayFastDir, ...
                                          CaImagingLongDelayFastFileList, ...
                                          params.minNumTrialToAnalysis, params); 
nonActiveNeuronIndex   = findNonActiveNeurons(nDataSet, params);
DataSetList(4).name    = 'Shuffle_Ca_Fast_Long_Delay';
DataSetList(4).params  = params; 
DataSetList(4).ActiveNeuronIndex = ~nonActiveNeuronIndex;
save([TempDatDir DataSetList(4).name '.mat'], 'nDataSet');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% long Ca slow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
nDataSet               = getCaImagingData(CaImagingLongDelaySlowDir, ...
                                          CaImagingLongDelaySlowFileList, ...
                                          params.minNumTrialToAnalysis, params); 
nonActiveNeuronIndex   = findNonActiveNeurons(nDataSet, params);
DataSetList(5).name    = 'Shuffle_Ca_Slow_Long_Delay';
DataSetList(5).params  = params; 
DataSetList(5).ActiveNeuronIndex = ~nonActiveNeuronIndex;
save([TempDatDir DataSetList(5).name '.mat'], 'nDataSet');

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
nDataSet               = getCaImagingData(CaImagingShortDelaySlowVirusDir, ...
                                          CaImagingShortDelaySlowVirusFileList, ...
                                          params.minNumTrialToAnalysis, params);
nonActiveNeuronIndex   = findNonActiveNeurons(nDataSet, params);
DataSetList(6).name    = 'Shuffle_Ca_Slow_Short_Delay_Virus';
DataSetList(6).params  = params; 
DataSetList(6).ActiveNeuronIndex = ~nonActiveNeuronIndex;
save([TempDatDir DataSetList(6).name '.mat'], 'nDataSet');


fileList            = {SpikeFileList; CaImagingShortDelayFastFileList; ...
                       CaImagingShortDelaySlowFileList;...
                       CaImagingLongDelayFastFileList; ...
                       CaImagingLongDelaySlowFileList; ...
                       CaImagingShortDelaySlowVirusFileList;...
                       SpikeFileList};


for nData           = 1:length(fileList)
    load([TempDatDir DataSetList(nData).name '.mat']);
    DataSetList(nData).cellinfo  = repmat(struct('fileName',1, 'nUnit', 1, ...
                                'AP_axis',1, 'ML_axis', 1, 'depth', 1,...
                                'expression','none', 'cellType',-1),length(nDataSet), 1);
    DataSetList(nData).ROCIndex  = ROCPop(nDataSet, DataSetList(nData).params);
    for nUnit  = 1:length(nDataSet)
        nFileList                                     = fileList{nData};
        DataSetList(nData).cellinfo(nUnit).fileName   = nFileList(nDataSet(nUnit).sessionIndex).name;
        DataSetList(nData).cellinfo(nUnit).nUnit      = nDataSet(nUnit).nUnit;
        DataSetList(nData).cellinfo(nUnit).AP_axis    = nDataSet(nUnit).AP_in_um;
        DataSetList(nData).cellinfo(nUnit).ML_axis    = nDataSet(nUnit).ML_in_um;
        DataSetList(nData).cellinfo(nUnit).depth      = nDataSet(nUnit).depth_in_um;
        DataSetList(nData).cellinfo(nUnit).expression = DataSetList(nData).params.expression;
        if strcmp(nDataSet(nUnit).cell_type, 'putative_interneuron')
            DataSetList(nData).cellinfo(nUnit).cellType   = 0;
        elseif strcmp(nDataSet(nUnit).cell_type, 'putative_pyramidal')
            DataSetList(nData).cellinfo(nUnit).cellType   = 1;
        end        
    end
end

save([TempDatDir 'DataListShuffle.mat'], 'DataSetList');
