%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Raw activity of all neurons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('../Func');
setDir;
load ([TempDatDir 'DataListC2SModel.mat']);

ylabels                 = {'dF/F', 'dF/F', 'dF/F', 'dR/R' };

if ~exist([PlotDir 'C2SModel'],'dir')
    mkdir([PlotDir 'C2SModel'])
end


for nData             = 1:length(DataSetList)
    load([TempDatDir DataSetList(nData).name '.mat'])
    plotMeanActivityImagescRasterOnly(nDataSet, DataSetList(nData).params, [], [], ylabels{nData}); 
    setPrint(6*2, 3*3, [PlotDir 'C2SModel/SingleUnitsImagescRasterOnly_' DataSetList(nData).name], 'svg')
end



close all;