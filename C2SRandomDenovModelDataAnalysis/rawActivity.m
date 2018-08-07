%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Raw activity of all neurons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('../Func');
setDir;
load ([TempDatDir 'DataListShuffle.mat']);
caDataSetList = DataSetList;
load ([TempDatDir 'DataListC2SRandomDeconvModel.mat']);

ylabels                 = {'dF/F', 'dF/F', 'dF/F', 'dR/R' };


for nData             = [1 2]
    load([TempDatDir caDataSetList(nData+2).name '_withOLRemoval.mat'])
    positivePeak = plotMeanActivityImagescRasterOnlyPositivePeak(nDataSet, DataSetList(nData).params, [], []); 
    load([TempDatDir DataSetList(nData).name '.mat'])
    C2SDataSet = nDataSet;
    C2SDataSet = C2SDataSet(~neuronRemoveList);
    C2SDataSet = C2SDataSet(positivePeak);
    plotMeanActivityImagescRasterOnlyC2S(C2SDataSet, DataSetList(nData).params, [], [], ylabels{nData}); 
    setPrint(6*2, 3*3, [PlotDir 'SingleUnitsImagescWithSort/SingleUnitsImagescRasterOnly_' DataSetList(nData).name '_withOLRemoval'])
end

load ([TempDatDir 'DataListShuffle.mat']);
caDataSetList = DataSetList;
load ([TempDatDir 'DataListC2SRandomDeconvModel.mat']);
nData = 10;
load([TempDatDir caDataSetList(nData).name '_withOLRemoval.mat'])
positivePeak = plotMeanActivityImagescRasterOnlyPositivePeak(nDataSet, caDataSetList(nData).params, [], []); 
load([TempDatDir DataSetList(3).name '.mat'])
C2SDataSet = nDataSet;
% C2SDataSet = C2SDataSet(~neuronRemoveList);
C2SDataSet = C2SDataSet(positivePeak);
plotMeanActivityImagescRasterOnlyC2S(C2SDataSet, DataSetList(3).params, [], [], 'DF/F'); 
setPrint(6*2, 3*3, [PlotDir 'SingleUnitsImagescWithSort/SingleUnitsImagescRasterOnly_' DataSetList(3).name '_withOLRemoval'])

close all;