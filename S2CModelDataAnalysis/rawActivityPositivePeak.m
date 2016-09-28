%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Raw activity of all neurons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('../Func');
setDir;
load ([TempDatDir 'DataListS2CModel.mat']);

load([TempDatDir 'Shuffle_Spikes.mat'])
positivePeak = plotMeanActivityImagescRasterOnlyPositivePeak (nDataSet, DataSetList(1).params, [], []);


for nData             = [3 4]
    load([TempDatDir DataSetList(nData).name '.mat'])
    plotMeanActivityImagescRasterOnly(nDataSet(positivePeak), DataSetList(nData).params, [], [], ''); 
    setPrint(6*2, 3*3, [PlotDir 'SingleUnitsImagescWithSort/SingleUnitsImagescRasterOnlyPositivePeak_' DataSetList(nData).name])
end



close all;