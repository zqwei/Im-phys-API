%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Raw activity of all neurons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('../Func');
setDir;
load ([TempDatDir 'DataListS2C6fModel.mat']);

load([TempDatDir 'Shuffle_Spikes_Nuo_Short_Delay.mat'])
depth       = [nDataSet.depth_in_um];
positivePeak = plotMeanActivityImagescRasterOnlyPositivePeak (nDataSet(depth<471), DataSetList(1).params, [], []);
setPrint(6*2, 3*3, [PlotDir 'SingleUnitsImagescWithSort/SingleUnitsImagescRasterMatched6f_' DataSetList(nData).name])

for nData             = 2
    load([TempDatDir DataSetList(nData).name '.mat'])
    plotMeanActivityImagescRasterOnly(nDataSet(positivePeak), DataSetList(nData).params, [], [], ''); 
    setPrint(6*2, 3*3, [PlotDir 'SingleUnitsImagescWithSort/SingleUnitsImagescRasterMatched6f_' DataSetList(nData).name])
end



% close all;