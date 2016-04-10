%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Raw activity of all neurons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('../Func');
setDir;
load ([TempDatDir 'DataListEphys.mat']);

ylabels                 = {'Firing rate', 'Firing rate', 'Firing rate', 'Firing rate', 'Firing rate', 'Firing rate'};

if ~exist([PlotDir 'EphysBinningTest'],'dir')
    mkdir([PlotDir 'EphysBinningTest'])
end


for nData             = 1:length(DataSetList)
    load([TempDatDir DataSetList(nData).name '.mat'])
    plotMeanActivityImagescRasterOnly(nDataSet, DataSetList(nData).params, [], [], ylabels{nData}); 
    setPrint(6*2, 3*3, [PlotDir 'EphysBinningTest/SingleUnitsImagescRasterOnly_' DataSetList(nData).name], 'svg')
end



close all;