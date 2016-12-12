%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Raw activity of all neurons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('../Func');
setDir;
load ([TempDatDir 'DataListShuffle.mat']);
load([TempDatDir 'Shuffle_Spikes.mat'])
positivePeak = plotMeanActivityImagescRasterOnlyPositivePeak (nDataSet, DataSetList(1).params, [], []);

load('DeconvGPResults/S2CC2SMCMCSingleTrial_FineTunedModeled_GP43.mat', 'nDataSet')

for nUnit        = 1:length(nDataSet)
    nDataSet(nUnit).unit_yes_trial      = nDataSet(nUnit).mcmc_yes_trial;
    nDataSet(nUnit).unit_no_trial       = nDataSet(nUnit).mcmc_no_trial;
end

nData = 1;
plotMeanActivityImagescRasterOnlyC2S(nDataSet(positivePeak), DataSetList(nData).params, [], [], ''); 
setPrint(6*2, 3*3, 'SingleUnitsImagescRasterOnlyPositivePeak')


close all;