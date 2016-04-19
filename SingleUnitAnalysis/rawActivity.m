%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Raw activity of all neurons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('../Func');
setDir;
load ([TempDatDir 'DataListShuffle.mat']);

ylabels                 = {'Fring Rate (Hz)', 'dF/F', 'dF/F', 'dF/F', 'dR/R' };
yAxes_set               = [0 60; -0.5 2.0; -0.5 2.0; -0.5 2.0; -0.5 2.0 ; -0.5 2.0];
lowFiringThres          = [15, 0.3, 0.3, 0.3, 0.3, 0.3];

if ~exist([PlotDir 'SingleUnitsImagescWithSort'],'dir')
    mkdir([PlotDir 'SingleUnitsImagescWithSort'])
end
% for nData             = 1:length(DataSetList)
%     load([TempDatDir DataSetList(nData).name '.mat'])
%     plotMeanActivityImagescWithSortWithCellinfo(nDataSet, DataSetList(nData).params, [], [], ylabels{nData}, lowFiringThres(nData), yAxes_set(nData,:)); 
%     setPrint(6*4, 3*3, [PlotDir 'SingleUnitsImagescWithSort/SingleUnitsImagescWithSort_' DataSetList(nData).name], 'tif')
% end
% 
% 
% for nData             = 1:length(DataSetList)
%     load([TempDatDir DataSetList(nData).name '.mat'])
%     plotMeanActivityImagescWithSortWithCellinfo(nDataSet(DataSetList(nData).ActiveNeuronIndex), DataSetList(nData).params, [], [], ylabels{nData}, lowFiringThres(nData), yAxes_set(nData,:)); 
%     setPrint(6*4, 3*3, [PlotDir 'SingleUnitsImagescWithSort/SingleActiveUnitsImagescWithSort_' DataSetList(nData).name], 'tif')
% end

for nData             = 1:length(DataSetList)-1
    load([TempDatDir DataSetList(nData).name '.mat'])
    plotMeanActivityImagescRasterOnly(nDataSet, DataSetList(nData).params, [], [], ylabels{nData}); 
    setPrint(6*2, 3*3, [PlotDir 'SingleUnitsImagescWithSort/SingleUnitsImagescRasterOnly_' DataSetList(nData).name])
end



close all;