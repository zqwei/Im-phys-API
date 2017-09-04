%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Raw activity of all neurons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('../Func');
setDir;
load ([TempDatDir 'DataListS2CGP43Model.mat']);

ylabels                 = {'dF/F', 'dF/F', 'dF/F', 'dF/F' };

if ~exist([PlotDir 'S2CGP43Model'],'dir')
    mkdir([PlotDir 'S2CGP43Model'])
end


for nData             = 2
    load([TempDatDir DataSetList(nData).name '.mat'])
    plotMeanActivityImagescRasterOnlyPositivePeak(nDataSet, DataSetList(nData).params, [], []); 
    setPrint(6*2, 3*3, [PlotDir 'S2CGP43Model/SingleUnitsImagescRasterOnlyPositivePeak_' DataSetList(nData).name])
end



close all;