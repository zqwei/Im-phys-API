%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 0.  Generate raw activity of all neurons sorted in different ways
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison_V2_0

addpath('../Func');
setDir;

minNumTrialToAnalysis  = 20;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % short Ca fast
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
% nDataSet               = getCaImagingDataNoNeuropil(CaImagingShortDelayFastDir, ...
%                                           CaImagingShortDelayFastFileList, ...
%                                           params.minNumTrialToAnalysis, params); 
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
nDataSet               = getCaImagingDataNoNeuropil(CaImagingShortDelaySlowDir, ...
                                          CaImagingShortDelaySlowFileList, ...
                                          params.minNumTrialToAnalysis, params); 
nonActiveNeuronIndex   = findNonActiveNeurons(nDataSet, params);
DataSetList(1).name    = 'NoNeuropil_Ca_Slow_Short_Delay';
DataSetList(1).params  = params; 
DataSetList(1).ActiveNeuronIndex = ~nonActiveNeuronIndex;
save([TempDatDir DataSetList(1).name '.mat'], 'nDataSet');

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
nDataSet               = getCaImagingDataNoNeuropil(CaImagingShortDelaySlowVirusDir, ...
                                          CaImagingShortDelaySlowVirusFileList, ...
                                          params.minNumTrialToAnalysis, params);
nonActiveNeuronIndex   = findNonActiveNeurons(nDataSet, params);
DataSetList(2).name    = 'NoNeuropil_Ca_Slow_Short_Delay_Virus';
DataSetList(2).params  = params; 
DataSetList(2).ActiveNeuronIndex = ~nonActiveNeuronIndex;
save([TempDatDir DataSetList(2).name '.mat'], 'nDataSet');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fileList            = {CaImagingShortDelaySlowFileList;...
                       CaImagingShortDelaySlowVirusFileList};


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

save([TempDatDir 'DataListNoNeuropil.mat'], 'DataSetList');
