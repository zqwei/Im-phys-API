%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Raw activity of all neurons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('../Func');
setDir;
load ([TempDatDir 'DataListS2CModel.mat']);
load([TempDatDir 'Shuffle_Spikes_Nuo_Short_Delay.mat'])
depth          = [nDataSet.depth_in_um];
spikeDataSet   = nDataSet(depth<471);
positivePeak = plotMeanActivityImagescRasterOnlyPositivePeak (spikeDataSet, DataSetList(1).params, [], []);

fName = 'Deconv_Ca_Fast_SShort_Delay';
load(['S2CC2S_' fName '.mat'])  
spkDataSet = nDataSet;
plotMeanActivityImagescRasterOnlyC2S(spkDataSet(positivePeak), DataSetList(1).params, [], [], ''); 
% plotMeanActivityImagescRasterOnlyC2S(spkDataSet, DataSetList(1).params, [], [], ''); 
setPrint(6*2, 3*3, [PlotDir 'SingleUnitsImagescWithSort/SingleUnitsImagescRasterOnlyPositivePeakS2CC2S_' fName])
