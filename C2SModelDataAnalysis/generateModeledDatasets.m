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
%
% Raw activity
%
% Spike
nDataSet               = getSpikeDataWithEphysTime(SpikingDataDir, SpikeFileList, params.minNumTrialToAnalysis, params.timeSeries, params.binsize);
spikeDataSet           = nDataSet;
% nDataSet               = getDFFSpike(spikeDataSet, params);
% nData                      = 1;
% DataSetList(nData).name    = 'Modeled_Spikes';
% DataSetList(nData).params  = params; 
ActiveNeuronIndex = findHighFiringUnits(spikeDataSet, params, minFiringRate);
% save([TempDatDir DataSetList(nData).name '.mat'], 'nDataSet');


% Modeled GCaMP6s
truncatedNormal        = truncate(makedist('Normal'), -1.5, 1.5);
params.Fm              = random(truncatedNormal, length(nDataSet), 1) * 10.6673 + 21.3240;
params.K               = random(truncatedNormal, length(nDataSet), 1) *  5.8337 + 13.9248;
params.n               = random(truncatedNormal, length(nDataSet), 1) *  0.4106 +  1.7531;
params.tau_r           = random(truncatedNormal, length(nDataSet), 1) *  0.0344 +  0.0728;
params.tau_d           = random(truncatedNormal, length(nDataSet), 1) *  0.5055 +  1.4551;
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

% Modeled GCaMP6s
truncatedNormal        = truncate(makedist('Normal'), -1.5, 1.5);
params.Fm              = random(truncatedNormal, length(nDataSet), 1) * 10.6673 + 21.3240;
params.K               = random(truncatedNormal, length(nDataSet), 1) *  5.8337 + 13.9248;
params.n               = random(truncatedNormal, length(nDataSet), 1) *  0.4106 +  1.7531;
params.tau_r           = random(truncatedNormal, length(nDataSet), 1) *  0.0344 +  0.0728;
params.tau_d           = random(truncatedNormal, length(nDataSet), 1) *  0.5055 +  1.4551;
params.Fm              = params.Fm;
params.tau_r           = params.tau_r;
params.tau_d           = params.tau_d/7;
params.intNoise        = 1.5;
params.extNoise        = 1.5;
nDataSet               = getFakeCaImagingData(spikeDataSet, params);
nData                      = 2;
DataSetList(nData).name    = 'Modeled_Ca_Short_Decay_No_Noise';
DataSetList(nData).params  = params; 
DataSetList(nData).ActiveNeuronIndex = ActiveNeuronIndex;
save([TempDatDir DataSetList(nData).name '.mat'], 'nDataSet'); 

save([TempDatDir 'DataListModeled.mat'], 'DataSetList');

for nData           = 1:length(DataSetList)
    load([TempDatDir DataSetList(nData).name '.mat']);
    DataSetList(nData).cellinfo  = repmat(struct('fileName',1, 'nUnit', 1, ...
                                'AP_axis',1, 'ML_axis', 1, 'depth', 1,...
                                'expression','none', 'cellType',-1),length(nDataSet), 1);
    for nUnit  = 1:length(nDataSet)
        nFileList                                     = SpikeFileList;
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

save([TempDatDir 'DataListC2SModel.mat'], 'DataSetList');

% 
% nData                       = 1;
% load([TempDatDir DataSetList(nData).name '.mat']);
% % only considering pyramidal cells
% selectedNeuronalIndex       = [DataSetList(nData).cellinfo(:).cellType]==1;
% minRate                     = 5; % Hz
% perMinRate                  = 0.4;
% ROCThres                    = 0.6;
% HighFRDataSet               = filterOutLowFR(nDataSet, DataSetList(nData).params, minRate, perMinRate, ROCThres);
% selectedNeuronalIndex       = selectedNeuronalIndex' & HighFRDataSet;
% save([TempDatDir 'DataListModeled.mat'], 'selectedNeuronalIndex', '-append');
