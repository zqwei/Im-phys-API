%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Single neuron choice probability index distribution (ROC), across
% trial periods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../Func');
setDir;
load ([TempDatDir 'DataListC2SModel.mat']);

if ~exist([PlotDir 'C2SModel'],'dir')
    mkdir([PlotDir 'C2SModel'])
end


for nData             = 1:length(DataSetList)
    load([TempDatDir DataSetList(nData).name '.mat'])
    plotROCPopAccLines(nDataSet, DataSetList(nData).params); 
    set(gca, 'TickDir', 'out')
    setPrint(8, 6, [PlotDir 'C2SModel/SingleUnitsROC_' DataSetList(nData).name], 'svg')
end

for nData             = 1:length(DataSetList)
    load([TempDatDir DataSetList(nData).name '.mat'])
    plotROCPopAccLines(nDataSet(DataSetList(nData).ActiveNeuronIndex), DataSetList(nData).params); 
    set(gca, 'TickDir', 'out')
    setPrint(8, 6, [PlotDir 'C2SModel/SingleActUnitsROC_' DataSetList(nData).name], 'svg')
end


close all