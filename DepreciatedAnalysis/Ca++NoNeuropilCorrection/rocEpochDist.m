%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Single neuron choice probability index distribution (ROC), across
% trial periods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../Func');
setDir;
load ([TempDatDir 'DataListNoNeuropil.mat']);

for nData             = 1:length(DataSetList)
    load([TempDatDir DataSetList(nData).name '.mat'])
    plotROCPopAccLines(nDataSet, DataSetList(nData).params); 
    set(gca, 'TickDir', 'out')
    setPrint(8, 6, [PlotDir 'NoNeuropil/SingleUnitsROC_' DataSetList(nData).name], 'svg')
end

for nData             = 1:length(DataSetList)
    load([TempDatDir DataSetList(nData).name '.mat'])
    plotROCPopAccLines(nDataSet(DataSetList(nData).ActiveNeuronIndex), DataSetList(nData).params); 
    set(gca, 'TickDir', 'out')
    setPrint(8, 6, [PlotDir 'NoNeuropil/SingleActUnitsROC_' DataSetList(nData).name], 'svg')
end

close all