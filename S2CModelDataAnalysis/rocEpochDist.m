%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Single neuron choice probability index distribution (ROC), across
% trial periods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../Func');
setDir;
load ([TempDatDir 'DataListS2CModel.mat']);

if ~exist([PlotDir 'S2CModel'],'dir')
    mkdir([PlotDir 'S2CModel'])
end


for nData             = 1:length(DataSetList)
    load([TempDatDir DataSetList(nData).name '.mat'])
    plotROCPopAccLines(nDataSet, DataSetList(nData).params); 
    set(gca, 'TickDir', 'out')
    setPrint(8, 6, [PlotDir 'S2CModel/SingleUnitsROC_' DataSetList(nData).name], 'svg')
end

for nData             = 1:length(DataSetList)
    load([TempDatDir DataSetList(nData).name '.mat'])
    plotROCPopAccLines(nDataSet(DataSetList(nData).ActiveNeuronIndex), DataSetList(nData).params); 
    set(gca, 'TickDir', 'out')
    setPrint(8, 6, [PlotDir 'S2CModel/SingleActUnitsROC_' DataSetList(nData).name], 'svg')
end


close all