%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Raw activity of all neurons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('../Func');
setDir;
load ([TempDatDir 'DataListModeled.mat']);

ylabels                 = {'Fring Rate (Hz)', 'dF/F', 'dF/F', 'dF/F', 'dR/R' };
yAxes_set               = [0 60; -0.5 2.0; -0.5 2.0; -0.5 2.0; -0.5 2.0 ; -0.5 2.0];
lowFiringThres          = [15, 0.3, 0.3, 0.3, 0.3, 0.3];

if ~exist([PlotDir 'ModeledSingleUnitsImagescWithSort'],'dir')
    mkdir([PlotDir 'ModeledSingleUnitsImagescWithSort'])
end
for nData             = 1:length(DataSetList)
    load([TempDatDir DataSetList(nData).name '.mat'])
    plotMeanActivityImagescWithSortWithCellinfo(nDataSet, DataSetList(nData).params, [], [], ylabels{nData}, lowFiringThres(nData), yAxes_set(nData,:)); 
    setPrint(6*4, 3*3, [PlotDir 'ModeledSingleUnitsImagescWithSort/SingleUnitsImagescWithSort_' DataSetList(nData).name], 'tif')
end


% for nData             = 1:length(DataSetList)
%     load([TempDatDir DataSetList(nData).name '.mat'])
%     plotMeanActivityImagescWithSortWithCellinfo(nDataSet(DataSetList(nData).ActiveNeuronIndex), DataSetList(nData).params, [], [], ylabels{nData}, lowFiringThres(nData), yAxes_set(nData,:)); 
%     setPrint(6*4, 3*3, [PlotDir 'ModeledSingleUnitsImagescWithSort/SingleActiveUnitsImagescWithSort_' DataSetList(nData).name], 'tif')
% end

close all;