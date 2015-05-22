%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1.  Raw activity of all neurons sorted in different ways
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison_V1_1
addpath('../Func');
setDir;
load ([TempDatDir 'DataListModeled.mat']);


ylabels                 = {'dR/R', 'dF/F', 'dF/F', 'dF/F', 'dF/F', 'dF/F'};
yAxes_set               = [0 60; -0.5 2.0; -0.5 2.0; -0.5 2.0; -0.5 2.0 ; -0.5 2.0];
lowFiringThres          = [15, 0.3, 0.3, 0.3, 0.3, 0.3];

if ~exist([PlotDir '/Single_UnitsImagescWithSortModeledData'],'dir')
    mkdir([PlotDir '/Single_UnitsImagescWithSortModeledData'])
end
% for nData             = 1:length(DataSetList)
%     load([TempDatDir DataSetList(nData).name '.mat'])
%     plotMeanActivityImagescWithSortWithCellinfo(nDataSet([DataSetList(nData).cellinfo(:).cellType]==1), DataSetList(nData).params, [], [], ylabels{nData}, lowFiringThres(nData), yAxes_set(nData,:)); 
%     setPrint(6*4, 3*3, [PlotDir '/Single_UnitsImagescWithSortModeledData/Single_UnitsImagescWithSort_' DataSetList(nData).name], 'tif')
% end

for nData             = 1:length(DataSetList)
    load([TempDatDir DataSetList(nData).name '.mat'])
    plotMeanActivityImagescWithSortWithCellinfo(nDataSet, DataSetList(nData).params, [], [], ylabels{nData}, lowFiringThres(nData), yAxes_set(nData,:)); 
    setPrint(6*4, 3*3, [PlotDir '/Single_UnitsImagescWithSortModeledData/Single_UnitsImagescWithSort_' DataSetList(nData).name], 'tif')
end

close all;

