%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Single neuron choice probability index distribution (ROC), across
% trial periods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../Func');
setDir;
load ([TempDatDir 'DataListShuffle.mat']);

if ~exist([PlotDir 'SingleUnitsROC'],'dir')
    mkdir([PlotDir 'SingleUnitsROC'])
end


for nData             = [3 4]
    load([TempDatDir DataSetList(nData).name '_withOLRemoval.mat'])
    plotROCPopLines(nDataSet(DataSetList(nData).ActiveNeuronIndex(~neuronRemoveList)), DataSetList(nData).params); 
    setPrint(8*3, 6, [PlotDir 'SingleUnitsROC/SingleActUnitsROCNoAcc_' DataSetList(nData).name '_withOLRemoval'])
end


close all