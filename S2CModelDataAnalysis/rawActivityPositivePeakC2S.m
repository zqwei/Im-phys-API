%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Raw activity of all neurons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('../Func');
setDir;
load ([TempDatDir 'DataListS2CModel.mat']);

load([TempDatDir 'Shuffle_Spikes.mat'])
positivePeak = plotMeanActivityImagescRasterOnlyPositivePeak (nDataSet, DataSetList(1).params, [], []);


for nData             = [3 4]
    load(['S2CC2S_' DataSetList(nData).name '.mat'])  
    spkDataSet = nDataSet;
    for nCell  = 1:length(nDataSet)
        spkDataSet(nCell).unit_yes_trial = spkDataSet(nCell).mcmc_yes_trial;
        spkDataSet(nCell).unit_no_trial  = spkDataSet(nCell).mcmc_no_trial;
    end
    plotMeanActivityImagescRasterOnly(spkDataSet(positivePeak), DataSetList(nData).params, [], [], ''); 
    setPrint(6*2, 3*3, [PlotDir 'SingleUnitsImagescWithSort/SingleUnitsImagescRasterOnlyPositivePeakS2CC2S_' DataSetList(nData).name])
end



close all;



load([TempDatDir 'Shuffle_Spikes_Nuo_Short_Delay.mat'])
depth          = [nDataSet.depth_in_um];
spikeDataSet   = nDataSet(depth<471);
positivePeak = plotMeanActivityImagescRasterOnlyPositivePeak (spikeDataSet, DataSetList(1).params, [], []);

fName = 'FineTunedModeled_GP517';
load(['S2CC2S_' fName '.mat'])  
spkDataSet = nDataSet;
for nCell  = 1:length(nDataSet)
    spkDataSet(nCell).unit_yes_trial = spkDataSet(nCell).mcmc_yes_trial;
    spkDataSet(nCell).unit_no_trial  = spkDataSet(nCell).mcmc_no_trial;
end
plotMeanActivityImagescRasterOnly(spkDataSet(positivePeak), DataSetList(1).params, [], [], ''); 
setPrint(6*2, 3*3, [PlotDir 'SingleUnitsImagescWithSort/SingleUnitsImagescRasterOnlyPositivePeakS2CC2S_' fName])
