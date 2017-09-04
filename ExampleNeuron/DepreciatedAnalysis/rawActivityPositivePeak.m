%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Raw activity of all neurons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('../Func');
setDir;
load ([TempDatDir 'DataListShuffle.mat']);
load([TempDatDir 'Shuffle_Spikes.mat'])
positivePeak = plotMeanActivityImagescRasterOnlyPositivePeak (nDataSet, DataSetList(1).params, [], []);

load('S2CC2SNL.mat', 'nDataSet')

for nUnit        = 1:length(nDataSet)
    nDataSet(nUnit).unit_yes_trial      = nDataSet(nUnit).mcmc_yes_trial(:, 56:132);
    nDataSet(nUnit).unit_no_trial       = nDataSet(nUnit).mcmc_no_trial(:, 56:132);
end

nData = 1;
plotMeanActivityImagescRasterOnly(nDataSet(positivePeak), DataSetList(nData).params, [], [], ''); 
setPrint(6*2, 3*3, 'SingleUnitsImagescRasterOnlyPositivePeak_MCMC')


close all;