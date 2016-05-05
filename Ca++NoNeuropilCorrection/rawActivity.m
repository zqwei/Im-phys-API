%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Raw activity of all neurons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('../Func');
setDir;
load ([TempDatDir 'DataListNoNeuropil.mat']);

ylabels                 = {'dF/F', 'dF/F'};
yAxes_set               = [-0.5 2.0; -0.5 2.0];
lowFiringThres          = [0.3, 0.3];

if ~exist([PlotDir 'NoNeuropil'],'dir')
    mkdir([PlotDir 'NoNeuropil'])
end
for nData             = 1:length(DataSetList)
    load([TempDatDir DataSetList(nData).name '.mat'])
    plotMeanActivityImagescRasterOnly(nDataSet, DataSetList(nData).params, [], [], ylabels{nData}); 
    setPrint(6*2, 3*3, [PlotDir 'NoNeuropil/SingleUnitsImagescRasterOnly_' DataSetList(nData).name], 'svg')
end
close all;