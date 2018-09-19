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

minNumTrialToAnalysis  = 15;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% short Ca slow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params.frameRate       =  29.68/2;
params.binsize         =  1/params.frameRate;
params.polein          =  -2.4;
params.poleout         =  -1.2;
start_time             =  params.polein - 0.5;
minTimeToAnalysis      =  round(start_time * params.frameRate);
maxTimeToAnalysis      =  round(1.2 * params.frameRate);
params.timeWindowIndexRange  = minTimeToAnalysis : maxTimeToAnalysis;
params.timeSeries      = params.timeWindowIndexRange * params.binsize;
params.minNumTrialToAnalysis =  minNumTrialToAnalysis;
params.expression      = 'Transgentic';

load([TempDatDir 'Shuffle_Ca_Fast_SShort_Delay_withOLRemoval.mat'])

nDataSet               = getFakeSpikeSingleTrialOOPSIData(nDataSet);  
nonActiveNeuronIndex   = findNonActiveNeurons(nDataSet, params);
nData                  = 1;
DataSetList(nData).name    = 'ModelSpikeMCMCSingleTrial_OOPSI_Ca_Fast_SShort_Delay';
DataSetList(nData).params  = params; 
DataSetList(nData).ActiveNeuronIndex = ~nonActiveNeuronIndex;
save([TempDatDir DataSetList(nData).name '.mat'], 'nDataSet');

save([TempDatDir 'DataListC2SMCMCSingleTrial6fModel.mat'], 'DataSetList');



load([TempDatDir 'FineTunedModeled_GP517.mat'])
for nCell = 1:length(nDataSet)
    nUnitData     = nDataSet(nCell).unit_yes_trial;
    nDataSet(nCell).mcmc_yes_trial = imagingToSpike(nUnitData)/binsize;
    
    nUnitData     = nDataSet(nCell).unit_no_trial;
    nDataSet(nCell).mcmc_no_trial  = imagingToSpike(nUnitData)/binsize;
    
end

save(['S2CC2S_' 'FineTunedModeled_GP517.mat'], 'nDataSet');
